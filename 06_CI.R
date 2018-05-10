##### Contains the code used to evaluate the performance of confidence intervals

rm(list=ls())
library(GenomicRanges); library(rtracklayer); library(stringr); library(Biostrings); library(Rsamtools); 
library(Biostrings); library(GenomicFeatures); library(BSgenome.Hsapiens.UCSC.hg38)
library(recountNNLS); library(recountNNLSdata)

##########################################################################################
### Import the annotation and subset to simulate from lengths > 150bp txs on chr1
##########################################################################################
setwd("/dcl01/leek/data/ta_poc/")
library(BSgenome.Hsapiens.UCSC.hg38)

### Select the transcripts to simulate
hg38 = BSgenome.Hsapiens.UCSC.hg38
TxDb = makeTxDbFromGFF("GencodeV25/gencodeV25.coding.chr1.gff3")
	TxL = transcriptLengths(TxDb)
	TxSeq = extractTranscriptSeqs(hg38, TxDb)
	names(TxSeq) = TxL$tx_name
	TxL = TxL[TxL$tx_len>=150,]
	TxSeq = TxSeq[match(TxL$tx_name, names(TxSeq))]

	set.seed(1337)
	genes = unique(TxL$gene_id)
	txs = sample(TxL$tx_name, 2000, replace=F)
	TxL_sim = TxL[TxL$tx_name %in% txs, ]
	TxSeq_sim = TxSeq[names(TxSeq) %in% txs, ]

	samp_size=100
	samps = paste0("sample_", c(paste0("0", 1:9), 10:samp_size))

### Perform the simulation studies for each condition
rls = c(37, 50, 75, 100, 150)
for(paired in 1:2){
	for(rl in rls){
		##########################################################################################
		### Simulate 20x coverage using polyester for transcripts
		##########################################################################################
		if(paired==2){
			outdir=paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/ci_tuning/", rl, "_2")
		}else{
			outdir=paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/ci_tuning/", rl)
		}
		dir.create(outdir); setwd(outdir)
		writeXStringSet(TxSeq_sim, "tx.fasta")

		reads = round(width(TxSeq_sim)*20/rl/paired)
		count_mat = matrix(rep(reads, samp_size), ncol=samp_size, byrow=F)
		rownames(count_mat) = TxL_sim$tx_name
		if(paired==2){
			save(count_mat, file=paste0("~/CI_cal_countmat_", rl, "_2.rda"))
		}else{
			save(count_mat, file=paste0("~/CI_cal_countmat_", rl, ".rda"))
		}
		polyester::simulate_experiment_countmat('tx.fasta', 
			readmat=count_mat, outdir="reads", paired=(paired==2), readlen=rl)	

		##########################################################################################
		### Align using rail-rna
		##########################################################################################
		if(paired==1){
			outdir=paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/ci_tuning/", rl)
			setwd(outdir)
			samples = paste0("sample_", c(paste0("0", 1:9), 10:100))
			man = cbind(paste0(outdir, "/reads/", samples, ".fasta"), 0, samples)
		}else{
			outdir=paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/ci_tuning/", rl, "_2/reads/")
			setwd(outdir)
			samples = paste0("sample_", c(paste0("0", 1:9), 10:100))
			man = cbind(paste0(outdir, samples, "_1.fasta"), 0, 
				paste0(outdir, samples, "_2.fasta"), 0, samples)
		}
		rail_dir = paste0(outdir, "/reads/rail")
		dir.create(rail_dir); setwd(rail_dir)
		write.table(man, sep="\t", quote=F, col.names=F, row.names=F, file="manifest.txt")	
		index="/dcl01/leek/data/HG38_rail/Homo_sapiens/UCSC/hg38/Sequence/"
		system2("rail-rna", paste0("go local --force --deliverables tsv bw -x ", index, "/BowtieIndex/genome ", index, "/Bowtie2Index/genome -m manifest.txt"))	

		##########################################################################################
		### Quantify with recountNNLS
		##########################################################################################
		if(paired==1){
			outdir=paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/ci_tuning/", rl)	
		}else{
			outdir=paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/ci_tuning/", rl, "_2")	
		}
		bw = paste0(outdir, '/reads/rail/rail-rna_out/coverage_bigwigs/', samps, '.bw')
		jx_file = paste0(outdir, '/reads/rail/rail-rna_out/cross_sample_results/junctions.tsv.gz')
		table = data.frame(project=rl, run=samps, bigwig_path=bw, rls=rl, paired_end=(paired==2))
		pheno = processPheno(table)
		rse_tx = recountNNLS(pheno, jx_file=jx_file, cores=20)
		if(paired==2){
			save(rse_tx, file=paste0("~/CI_", rl, "_2.rda"))	
		}else{
			save(rse_tx, file=paste0("~/CI_", rl, ".rda"))	
		}
	}
}

#### Code for figure showing number of transcripts that attains 
### 95% coverage with different read lengths
rm(list=ls())
means_list = list()
for(paired in 1:2){
	for(i in 1:length(rls)){
		rl = rls[i]
		if(paired==1){
			load(paste0("~/CI_cal_countmat_", rl, ".rda"))
			load(paste0("~/CI_", rl, ".rda"))	
		}else{
			load(paste0("~/CI_cal_countmat_", rl, "_2.rda"))
			load(paste0("~/CI_", rl, "_2.rda"))	
		}
		rse_sub = rse_tx[match(rownames(count_mat), rownames(rse_tx)),]
		bs = assays(rse_sub)$fragments
		se = assays(rse_sub)$ses
		score = assays(rse_sub)$scores
		df = assays(rse_sub)$df
		stats = qt(0.975, as.numeric(df[,1]))
		cil = bs-stats*se
		ciu = bs+stats*se
		hit = (cil<=count_mat)*(ciu>=count_mat)

		hits = apply(hit, 1, mean)
		score_cat = cut(score[,1], seq(0, 1, by=0.1))

		thresholds = seq(0.9, 1, by=0.01)
		means = NULL
		for(thresh in thresholds){
			means = rbind(means,  by(hits, score_cat, function(x) mean(x<thresh, na.rm=T)))
		}
		rownames(means) = thresholds
		means_list[[(i+(paired-1)*5)]] = means
	}
}

### Ploting
png(paste0('~/graphics/CI_cal.png'), width=900, height=300, type="cairo")
par(mfrow=c(1,3))
	# Panel A
	par(ps = 12, cex = 1, cex.main = 1)
	par(mar=c(4, 4, 4, 2))
	hist(score, breaks=seq(0, 1, by=0.1), xlab="Uniqueness", 
		main="Distribution of Uniqueness Scores", col=0)
	u = par()$usr; rect(u[1], u[3], u[2], u[4], col=rgb(0, 0, 0, 0.2), border=NA)
	text(x=u[2]-(u[2]-u[1])/20, y=u[4]-(u[4]-u[3])/20, label="(A)")
	hist(score, breaks=seq(0, 1, by=0.1), add=T, col=0)
	par(mar=c(4, 6, 4, 1))
	par(mar=c(4, 4, 4, 0))

	# Panel B
	ind=6
	plot(means_list[[1]][ind,], ylab="Proportion of txs w/o nominal coverage", 
		xlab="", xaxt="n", main="95% CI Single-End", pch=19, ty="b", ylim=c(0, max(means_list[[1]][6,])))
	for(i in 2:length(rls)){
		points(means_list[[i]][ind,], col=i, pch=19, ty="b")
	}
	u = par()$usr; rect(u[1], u[3], u[2], u[4], col=rgb(0, 0, 0, 0.2), border=NA)
	axis(side=1, at = 1:length(colnames(means)), labels = colnames(means), las=2)
	text(x=mean(u[1:2]), y=u[3], labels="Uniqueness", pos=3)
	text(x=u[1]+(u[2]-u[1])/20, y=u[3]+(u[4]-u[3])/20, label="(B)")
	legend('topright', legend = rls[1:3], col=1:3, fill=1:3, cex=0.8, ncol=3)

	# Panel C
	par(mar=c(4, 0, 4, 4))
	plot(means_list[[6]][ind,], ylab="", 
		xlab="", xaxt="n", main="95% CI Paired-End", pch=19, ty="b", ylim=c(0, max(means_list[[1]][6,])), yaxt="n")
	for(i in 2:length(rls)){
		points(means_list[[(i+5)]][ind,], col=i, pch=19, ty="b")
	}
	u = par()$usr; rect(u[1], u[3], u[2], u[4], col=rgb(0, 0, 0, 0.2), border=NA)
	axis(side=1, at = 1:length(colnames(means)), labels = colnames(means), las=2)
	axis(side=4)
	text(x=mean(u[1:2]), y=u[3], labels="Uniqueness", pos=3)
	text(x=u[1]+(u[2]-u[1])/20, y=u[3]+(u[4]-u[3])/20, label="(C)")
	legend('topleft', legend = rls[4:5], col=4:5, fill=4:5, ncol=2, cex=0.8)
dev.off()
