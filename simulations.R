##########################################################################################
module load cufflinks/2.2.1
module load python/2.6.6
module load bowtie/1.1.1
PATH=$PATH:/users/jmfu/kallisto_linux-v0.43.0:/users/jmfu/RSEM-1.3.0
##########################################################################################


##########################################################################################
### Preparing simulation truth (R)
##########################################################################################
rm(list=ls())
library(polyester); library(rtracklayer); library(Biostrings); library(stringr)
library(Biostrings); library(exomeCopy); library(GenomicFeatures); library(BSgenome.Hsapiens.UCSC.hg38)

### Obtain transcript sequences
hg38 = BSgenome.Hsapiens.UCSC.hg38
TxDb = makeTxDbFromGFF("/dcl01/leek/data/ta_poc/GencodeV25/gencodeV25.coding.gff3")
	TxL = transcriptLengths(TxDb)
	TxSeq = extractTranscriptSeqs(hg38, TxDb)
	names(TxSeq) = TxL$tx_name
	writeXStringSet(TxSeq, "/dcl01/leek/data/ta_poc/geuvadis/simulation/tx_full.fasta")

### Generate number of reads to simulate
library(gtools)
simp = function(gene, sim_info){
	message(which(genes==gene))
	sim_info_sub = sim_info[sim_info$gene_id == gene,,drop=FALSE]
	k = dim(sim_info_sub)[1]
	if(k>1){
		abund = rdirichlet(1, rep(1/k, k))
	}else{
		abund=1
	}	
	return(abund)
}
getTxs = function(gene, sim_info){
	message(which(genes==gene))
	sim_info_sub = sim_info[sim_info$gene_id == gene,,drop=FALSE]
	return(sim_info_sub$tx_name)
}
seed=76
	
	set.seed(seed)
	### Draw the relative expression of isoforms
	genes = unique(TxL$gene_id)
	tx_ps = mclapply(genes, simp, TxL, mc.cores=20)
	tx_ids = mclapply(genes, getTxs, TxL, mc.cores=20)
	tx_p = do.call(c, tx_ps)
	tx_id = do.call(c, tx_ids)

	set.seed(seed+1)
	### Draw the number reads for each gene
	num_tx = sapply(tx_ps, length)
	geneReads = rep(rnbinom(length(genes), 4, 0.01), times=num_tx)
	geneExpressed = rep(rbinom(length(genes), 1, 0.82), times=num_tx)
	geneNames = rep(genes, times=num_tx)
	
	### Combine gene and isoform statistics to determine number of reads per tx
	sim_info = data.frame(gene_id = geneNames, transcript_id = tx_id, geneExpressed = geneExpressed, geneReads=geneReads, tx_p = tx_p)
	sim_info = sim_info[match(TxL$tx_name, sim_info$transcript_id),]
	sim_info$transcriptsReads = round(sim_info$geneReads * sim_info$geneExpressed * sim_info$tx_p)
		test1 = sum(sim_info$transcript_id==names(TxSeq))
		test2 = sum(by(sim_info$tx_p, sim_info$gene_id, sum)>0.99999)
	sim_info$nexon = TxL$nexon
	sim_info$tx_len = TxL$tx_len

write.table(sim_info, row.names=F, file="/dcl01/leek/data/ta_poc/geuvadis/simulation/sim_info_full.txt")

### Simulate
exp_ind = sim_info$transcriptsReads>0
sim_info_expressed = sim_info[exp_ind,]
write.table(sim_info_expressed, row.names=F, file="/dcl01/leek/data/ta_poc/geuvadis/simulation/sim_info_full_expressed.txt")
TxSeq_expressed = TxSeq[exp_ind]


setwd('/dcl01/leek/data/ta_poc/geuvadis/simulation/')
conditions = c('37_1_full', '37_2_full', '50_1_full', '50_2_full', '75_2_full', '100_1_full', '100_2_full', '150_1_full', '150_2_full')
for(condition in conditions){
	dir.create(condition); setwd(condition)
	writeXStringSet(TxSeq_expressed, paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/", condition, "/tx_full_expressed.fasta"))
	dir.create("reads")
	readmat = matrix(sim_info_expressed$transcriptsReads, ncol=1)
	set.seed(seed+2)
	info = str_split_fixed(condition, "_", n=3)
	simulate_experiment_countmat('tx_full_expressed.fasta', readmat=readmat, outdir="reads", paired=(info[2]==2), readlen=as.numeric(info[1]))
	setwd('/dcl01/leek/data/ta_poc/geuvadis/simulation/')
}


##########################################################################################
### Run methods
##########################################################################################

	rl = 100
	paired = 1
	condition = paste0(rl, "_", paired, "_full")
	samples = paste0("sample_0", 1)
	outdir = paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/", condition)
	dir.create(outdir)
	setwd(outdir)

	##########################################################################################
	### Rail-RNA
	##########################################################################################
	if(paired==1){
		files = paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/", condition, "/reads/", samples, ".fasta")
		man = cbind(files, 0, samples)
		rail_dir = paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/", condition, "/rail")
		dir.create(rail_dir)
		write.table(man, sep="\t", quote=F, col.names=F, row.names=F, file=paste0(rail_dir, "/manifest.txt"))	
	}else{
		files = paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/", condition, "/reads/", samples)
		man = cbind(paste0(files, "_1.fasta"), 0, paste0(files, "_2.fasta"), 0, samples)
		rail_dir = paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/", condition, "/rail")
		dir.create(rail_dir)
		write.table(man, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, file=paste0(rail_dir, "/manifest.txt"))	
	}
	dir.create(paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/", condition, "/recount/"))
	dir.create(paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/", condition, "/recount/logs"))
	index="/dcl01/leek/data/HG38_rail/Homo_sapiens/UCSC/hg38/Sequence/"
	setwd(rail_dir)
	system2("rail-rna", paste0("go local -x ", index, "/BowtieIndex/genome ", index, "/Bowtie2Index/genome -m ", rail_dir, "/manifest.txt"))

	bw = paste0('/dcl01/leek/data/ta_poc/geuvadis/simulation/', condition, '/rail/rail-rna_out/coverage_bigwigs/sample_01.bw')
	jx_file = paste0('/dcl01/leek/data/ta_poc/geuvadis/simulation/', condition, '/rail/rail-rna_out/cross_sample_results')
	table = data.frame(project=rl, run=paired, bigwig_path=bw, rls=rl, paired_end=FALSE)
	pheno = processPheno(table)
	rse_tx = recountNNLS(pheno, jx_file=jx_file, cores=20)

	save(rse_tx, file=paste0('/dcl01/leek/data/ta_poc/geuvadis/simulation/', condition, '/recount/rse_tx.rda'))

	##########################################################################################
	### HISAT2-StringTie (2.05-1.3.3b)
	##########################################################################################
	setwd(outdir)
	dir.create("hisat2")
	dir.create("hisat2/aligned")
	dir.create("hisat2/sorted")
	for(s in samples){
		print(s); flush.console()
		outPath = paste0(outdir, "/hisat2/sorted/", s, ".ss.bam")
		readPath = paste0(outdir, "/reads/", s, ".fasta")
		outBam = paste0(outdir, "/hisat2/sorted/", s, ".ss")
		outSam = paste0(outdir, "/hisat2/aligned/", s, ".ss.sam")
		intermediateBam = paste0(outdir, "/hisat2/aligned/", s, "ss..bam")
		if(paired==1){
			system2("hisat2", paste0("--dta -f -x /dcl01/leek/data/ta_poc/hisat2-full-index/hisat2-ind-full -U ", readPath, " -S ", outSam))
		}else{
			system2("hisat2", paste0("--dta -f -x /dcl01/leek/data/ta_poc/hisat2-full-index/hisat2-ind-detailed -1 ", readPath, "_1.fasta", " -2 ", readPath, "_2.fasta ", "-S ", outSam))
		}
		system2("samtools", paste0("view -b ", outSam, " > ", intermediateBam))
		system2("samtools", paste0("sort ", intermediateBam, " ", outBam))
	}

	stringtie_dir = paste0(outdir, "/stringtie")
	annotation="/dcl01/leek/data/ta_poc/GencodeV25/gencodeV25.coding.gtf"
	system2("mkdir", stringtie_dir)
	system2("cd", stringtie_dir)

	for(s in samples){
		workingDir = paste0(stringtie_dir, "/", s)
		system2("mkdir", workingDir)
		system2("stringtie", paste0("-o ", workingDir, "/cov-e.gtf -G ", annotation, " -e /dcl01/leek/data/ta_poc/geuvadis/simulation/", condition, "/hisat2/sorted/", s, ".ss.bam") )
		system2("stringtie", paste0("-o ", workingDir, "/cov.gtf -G ", annotation, " /dcl01/leek/data/ta_poc/geuvadis/simulation/", condition, "/hisat2/sorted/", s, ".ss.bam") )
	}

	##########################################################################################
	### HISAT2-cufflinks (2.0.5)
	##########################################################################################
	setwd(outdir)
	dir.create("hisat2")
	dir.create("hisat2/aligned")
	dir.create("hisat2/sorted")

	for(s in samples){
		print(s); flush.console()
		outPath = paste0(outdir, "/hisat2/sorted/", s, ".cl.bam")
		readPath = paste0(outdir, "/reads/", s, ".fasta")
		outBam = paste0(outdir, "/hisat2/sorted/", s, ".cl")
		outSam = paste0(outdir, "/hisat2/aligned/", s, ".cl.sam")
		intermediateBam = paste0(outdir, "/hisat2/aligned/", s, "cl..bam")
		if(paired==1){
			system2("hisat2", paste0("--dta-cufflinks -f -x /dcl01/leek/data/ta_poc/hisat2-full-index/hisat2-ind-full -U ", readPath, " -S ", outSam))	
		}else{
			system2("hisat2", paste0("--dta-cufflinks -f -x /dcl01/leek/data/ta_poc/hisat2-full-index/hisat2-ind-detailed -1 ", readPath, "_1.fasta", " -2 ", readPath, "_2.fasta ", "-S ", outSam))
		}
		system2("samtools", paste0("view -b ", outSam, " > ", intermediateBam))
		system2("samtools", paste0("sort ", intermediateBam, " ", outBam))
	}

	cufflinks_dir = paste0(outdir, "/hisat2/cufflinks")
	system2("mkdir", cufflinks_dir)
	annotation="/dcl01/leek/data/ta_poc/GencodeV25/gencodeV25.coding.gtf"

	for(s in samples){
		workingDir = paste0(cufflinks_dir, "/", s)
		system2("mkdir", workingDir)
		system2("cufflinks", paste0("-m 250 -s 25 --total-hits-norm --no-effective-length-correction --no-length-correction -o ",
			workingDir, " -G ", annotation, " /dcl01/leek/data/ta_poc/geuvadis/simulation/", condition, "/hisat2/sorted/", s, ".cl.bam"))
	}

	##########################################################################################
	### Kallisto (0.43.0)
	##########################################################################################
	kallisto_dir = paste0(outdir, "/kallisto")
	system2("mkdir", kallisto_dir)
	index="/dcl01/leek/data/ta_poc/geuvadis/kallisto/kallisto_index"

	for(s in samples){
		workingDir = paste0(kallisto_dir, "/", s)
		system2("mkdir", workingDir)
		if(paired==1){
			system2("kallisto", paste0("quant -i ", index, " -o ", workingDir, " --single -l 250 -s 25 ", outdir, "/reads/", s, ".fasta"))
		}else{
			system2("kallisto", paste0("quant -i ", index, " -o ", workingDir , " ", outdir, "/reads/", s, "_1.fasta ", outdir, "/reads/", s, "_2.fasta"))
		}
	}

	##########################################################################################
	### Salmon (0.8.2)
	##########################################################################################
	index="/dcl01/leek/data/ta_poc/geuvadis/GencodeV25/salmon_index "
	salmon_dir = paste0(outdir, "/salmon")
	dir.create(salmon_dir)
	for(s in samples){
		if(paired==1){
			system2("/users/jmfu/Salmon-0.8.2_linux_x86_64/bin/salmon", 
				paste0("quant -i ", index, " -l A -r ", outdir, "/reads/", s, ".fasta -o ", salmon_dir, "/", s))	
		}
		else{
			system2("/users/jmfu/Salmon-0.8.2_linux_x86_64/bin/salmon", 
				paste0("quant -i ", index, " -l A -1 ", outdir, "/reads/", s, "_1.fasta -2 ", outdir, "/reads/", s, "_2.fasta -o ", salmon_dir, "/", s))	
		}
	}

	##########################################################################################
	### RSEM
	##########################################################################################
	rsem_dir = paste0(outdir, "/rsem")
	dir.create(rsem_dir)
	index = "/dcl01/leek/data/ta_poc/rsem-index/rsem "

	for(s in samples){
		if(paired==1){
			system2("rsem-calculate-expression", paste0("--no-qualities ", outdir, "/reads/", s, ".fasta ", index, rsem_dir, "/", s))
		}else{
			system2("rsem-calculate-expression", paste0("--no-qualities --paired-end ",  outdir, "/reads/", s, "_1.fasta ", outdir, "/reads/", s, "_2.fasta ", index, rsem_dir, "/", s))
		}
	}

##########################################################################################
### Compile results
##########################################################################################
	library(rtracklayer)
	truth = read.table("/dcl01/leek/data/ta_poc/geuvadis/simulation/sim_info_full2.txt", header=T)

	# load(paste0('/dcl01/leek/data/ta_poc/geuvadis/simulation/', condition, '/recount/rse_tx.rda'))
	recountNNLS = assays(rse_tx)$counts
	recountNNLS = recountNNLS[match(truth$transcript_id, rownames(recountNNLS)),]

	kallisto = read.table(paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/", condition, "/kallisto/sample_01/abundance.tsv"), header=T)
	kl_mat = match(truth$transcript_id, kallisto$target_id)
	kl = kallisto$est_counts[kl_mat]

	setwd(paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/", condition, "/stringtie"))
	info = import("sample_01/cov-e.gtf")
	info = info[info$type=="transcript"]
	st_cov = data.frame(transcript_id = info$transcript_id, info$cov)
	st_cov = st_cov[match(truth$transcript_id, st_cov$transcript_id),]
	st = as.numeric(as.character(st_cov[,2]))*truth$tx_len/75

	setwd(paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/", condition, "/salmon"))
	tab = read.table("sample_01/quant.sf", header=T)
	salmon = data.frame(transcript_id = tab[,1], reads = tab[,5])
	salmon = salmon[match(truth$transcript_id, salmon$transcript_id),]
	sl = salmon[,2]

	setwd(paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/", condition, "/rsem"))
	tab = read.table("sample_01.isoforms.results", header=T)
	rsem = cbind(tab[,5])
	rownames(rsem) = tab$transcript_id
	rsem = rsem[match(truth$transcript_id, rownames(rsem)),]

	tmp = import(paste0('/dcl01/leek/data/ta_poc/geuvadis/simulation/', condition, '/hisat2/cufflinks/sample_01/transcripts.gtf'))
	cufflinks = tmp[tmp$type=="transcript"]
	cufflinks_cov = data.frame(transcript_id = cufflinks$transcript_id, cufflinks$cov)
	cufflinks_cov = cufflinks_cov[match(truth$transcript_id, cufflinks_cov$transcript_id),]
	cl = as.numeric(as.character(cufflinks_cov[,2]))*truth$tx_len/75

	info = cbind(recountNNLS, st, kl, sl, rsem, cl)
	info[is.na(info)] = 0
	l2info = log(info+1, 2)
	l2truth = log(truth$transcriptsReads+1, 2)
	info_gene = matrix(0, ncol=6, nrow=19950)
	for(i in 1:6){
		info_gene[,i] = by(info[,i], truth$gene_id, sum)
	}
	truth_gene = as.numeric(by(truth$transcriptsReads, truth$gene_id, sum))
	l2info_gene = log(info_gene+1, 2)
	l2truth_gene = log(truth_gene+1, 2)
	
	save(info, truth, l2info, l2truth, file="~/sim_test.rda")
