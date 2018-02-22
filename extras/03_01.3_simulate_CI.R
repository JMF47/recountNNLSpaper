rm(list=ls())
library(GenomicRanges); library(rtracklayer); library(stringr); library(Biostrings); library(Rsamtools); 
library(Biostrings); library(GenomicFeatures); library(BSgenome.Hsapiens.UCSC.hg38)

##########################################################################################
### Import the annotation and subset to simulate from lengths > 150bp
##########################################################################################

setwd("/dcl01/leek/data/ta_poc/")
library(BSgenome.Hsapiens.UCSC.hg38)
hg38 = BSgenome.Hsapiens.UCSC.hg38
TxDb = makeTxDbFromGFF("GencodeV25/gencodeV25.coding.gff3")
	TxL = transcriptLengths(TxDb)
	TxSeq = extractTranscriptSeqs(hg38, TxDb)
		names(TxSeq) = TxL$tx_name
		TxL_sub = TxL[TxL$gene_id=="ENSG00000197147.12",]
		TxSeq_sub = TxSeq[match(TxL_sub$tx_name, names(TxSeq))]

##########################################################################################
### Simulate the reads using polyester
##########################################################################################
TxSeq = import('/dcl01/leek/data/ta_poc/geuvadis/simulation/rsem_based/100_1/tx.fasta')
TxL = rowData(rse_tx)
	TxL_sub = TxL[TxL$gene_id=="ENSG00000188257.10",]
	TxSeq_sub = TxSeq[match(TxL_sub$tx_name, names(TxSeq))]
rl=150
paired=FALSE
# gene = 'ENSG00000197147.12'
gene = 'ENSG00000188257.10'
outdir=paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/ci_tuning/", rl)
dir.create(outdir); setwd(outdir)
writeXStringSet(TxSeq_sub, paste0(gene, ".fasta"))
samp_size = 100
# count_mat = matrix(rep(c(0, 500, 0, 50), samp_size), ncol=samp_size, byrow=F)
count_mat = matrix(rep(c(50, 850, 0, 0, 0, 0, 50), samp_size), ncol=samp_size, byrow=F)
set.seed(rl+paired)
polyester::simulate_experiment_countmat(paste0(gene, '.fasta'), 
	readmat=count_mat, outdir=gene, paired=paired, readlen=rl)

##########################################################################################
### Quantify
##########################################################################################
dir = paste0(outdir, "/", gene)
setwd(dir)
samples = paste0("sample_", c(paste0("0", 1:9), 10:samp_size))
### Rail-RNA
files = paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/ci_tuning/", rl, "/", gene, "/", samples, ".fasta")
man = cbind(files, 0, samples)
rail_dir = paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/ci_tuning/", rl, "/", gene, "/rail")
dir.create(rail_dir)
write.table(man, sep="\t", quote=F, col.names=F, row.names=F, file=paste0(rail_dir, "/manifest.txt"))	
index="/dcl01/leek/data/HG38_rail/Homo_sapiens/UCSC/hg38/Sequence/"
system2("rail-rna", paste0("go local --force --deliverables tsv bw -x ", index, "/BowtieIndex/genome ", index, "/Bowtie2Index/genome -m ", rail_dir, "/manifest.txt"))

### recountNNLS
library(recountNNLS); library(recountNNLSdata)
samp_size=100
samples = paste0("sample_", c(paste0("0", 1:9), 10:samp_size))
bw = paste0(rail_dir, '-rna_out/coverage_bigwigs/', samples, '.bw')
jx_file = paste0(rail_dir, '-rna_out/cross_sample_results/junctions.tsv.gz')
table = data.frame(project=rl, run=samples, bigwig_path=bw, rls=rl, paired_end=FALSE)
pheno = processPheno(table)
power=0
cores = 20

message(Sys.time(), paste0(" ### Processing read length group: ", rl))
pheno = pheno[pheno$rls_group==rl,,drop=FALSE]
## Load emission probability matrices
message(Sys.time(), " # Setting up model covariates")
data(list=paste0("matrix_", rl), package = "recountNNLSdata")
matrix_list <- eval(parse(text=paste0("matrix_", rl)))
genes = names(matrix_list)
## Load feature counts
message(Sys.time(), " # Compiling feature counts")
counts = getCounts(pheno[pheno$rls_group==rl,], jx_file, cores = cores)

message(Sys.time(), " # Executing model")
load(paste0("~/", rl, "_g2l.rda"))
# locus = 16472
locus = 15913
loci = unique(g2l$locus)
outcome = list()
for(i in 1:samp_size){
	message(i)
	outcome[[i]] = inferReads(locus, matrix_list, counts[,i,drop=F], power=0)
	# outcome[[i]] = inferReads0("ENSG00000197147.12", list_1p, matrix_list, counts[,i,drop=F], power=0)
}
bs = sapply(outcome, function(x) x[[1]])
ses = sqrt(sapply(outcome, function(x) x[[2]]))
scores1 = sapply(outcome, function(x) x[[4]])
scores2 = sapply(outcome, function(x) x[[5]])


# truth = matrix(rep(c(0, 500, 50), samp_size), ncol=samp_size, byrow=F)
truth = matrix(rep(c(50, 850, 0, 0, 0, 0, 50), samp_size), ncol=samp_size, byrow=F)
cil = bs-1.96*ses
ciu = bs+1.96*ses
results = (cil<=truth & ciu>=truth)
apply(results, 1, mean)

	### Similarity score
P = matrix_150[[1]]

scores = sapply(1:dim(P)[2], function(i, P) {
	rsums = apply(P, 2, sum)
	min(apply(abs(P[,-i] - P[,i]), 2, sum)/(rsums[i]+rsums[-i]))}, P)


####
truth$scores = scores
vw = data.frame(round(truth$scores*100), truth$reads, round(truth$us), truth$gene_id)
	colnames(vw)= c("score", "reads", "truth", "gene")
	vw_bad = vw[vw[,1]<5,]
bw_truth = by(vw_bad[,2], vw_bad[,4], sum)
bw_us = by(vw_bad[,3], vw_bad[,4], sum)
collapse_bad = cbind(bw_truth, bw_us)
	collapse_bad = collapse_bad[complete.cases(collapse_bad),]

mrd_bad = abs(vw_bad[,2]-vw_bad[,3])/(vw_bad[,2]+vw_bad[,3])
	mean(mrd_bad, na.rm=T)

vw_good = vw[vw[,1]>5,]
mrd_good = abs(vw_good[,2] - vw_good[,3])/(vw_good[,2]+vw_good[,3])	
	mean(mrd_good, na.rm=T)





