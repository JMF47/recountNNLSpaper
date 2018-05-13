##### Contains the code used for confidence interval performance evaluation

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
			save(count_mat, 
				file=paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/CI/CI_cal_countmat_", rl, "_2.rda"))
		}else{
			save(count_mat, file=paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/CI/CI_cal_countmat_", rl, ".rda"))
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
			save(rse_tx, file=paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/CI/CI_", rl, "_2.rda"))	
		}else{
			save(rse_tx, file=paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/CI/CI_", rl, ".rda"))	
		}
	}
}


