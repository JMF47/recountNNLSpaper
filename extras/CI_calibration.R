rm(list=ls())
library(GenomicRanges); library(rtracklayer); library(stringr); library(Biostrings); library(Rsamtools); 
library(Biostrings); library(GenomicFeatures); library(BSgenome.Hsapiens.UCSC.hg38)

##########################################################################################
### Import the annotation and subset to simulate from lengths > 150bp
##########################################################################################

setwd("/dcl01/leek/data/ta_poc/")
library(BSgenome.Hsapiens.UCSC.hg38)
hg38 = BSgenome.Hsapiens.UCSC.hg38
TxDb = makeTxDbFromGFF("GencodeV25/gencodeV25.coding.chr1.gff3")
	TxL = transcriptLengths(TxDb)
	TxSeq = extractTranscriptSeqs(hg38, TxDb)
	names(TxSeq) = TxL$tx_name
	TxL = TxL[TxL$tx_len>=150,]
	TxSeq = TxSeq[match(TxL$tx_name, names(TxSeq))]

##########################################################################################
### Pick transcripts and simulate 20/5x coverage using polyester
##########################################################################################
### Sample 1 transcript from each gene
set.seed(1337)
genes = unique(TxL$gene_id)
txs = NULL
for(g in genes){
	TxL_sub = TxL[TxL$gene_id==g,]
	txs = c(txs, sample(TxL$tx_name, 1))
}
# save(txs, file="~/selectex_txs.rda")

TxL_sim = TxL[TxL$tx_name %in% txs, ]
TxSeq_sim = TxSeq[names(TxSeq) %in% txs, ]

outdir=paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/ci_tuning/r2/")
dir.create(outdir); setwd(outdir)
writeXStringSet(TxSeq_sim, "tx.fasta")

samp_size = 100
reads = round(width(TxSeq_sim)*40/150)
# reads = round(width(TxSeq_sim)*5/150)

count_mat = matrix(rep(reads, samp_size), ncol=samp_size, byrow=F)
rownames(count_mat) = TxL_sim$tx_name
save(count_mat, file="~/CI_cal_countmat.rda")

set.seed(1337)
polyester::simulate_experiment_countmat('tx.fasta', readmat=count_mat, outdir="reads40", paired=FALSE, readlen=150)
# polyester::simulate_experiment_countmat('tx.fasta', readmat=count_mat, outdir="reads5", paired=FALSE, readlen=150)

##########################################################################################
### Align using rail-rna
##########################################################################################
setwd(outdir)
samples = paste0("sample_", c(paste0("0", 1:9), 10:samp_size))
### Rail-RNA
files = paste0(outdir, "reads40/", samples, ".fasta")
# files = paste0(outdir, "reads5/", samples, ".fasta")
man = cbind(files, 0, samples)
rail_dir = paste0(outdir, "/reads40/rail")
# rail_dir = paste0(outdir, "/reads5/rail")
dir.create(rail_dir)
setwd(rail_dir)
write.table(man, sep="\t", quote=F, col.names=F, row.names=F, file="manifest.txt")	
index="/dcl01/leek/data/HG38_rail/Homo_sapiens/UCSC/hg38/Sequence/"
system2("rail-rna", paste0("go local --force -p 20 --deliverables tsv bw -x ", index, "/BowtieIndex/genome ", index, "/Bowtie2Index/genome -m manifest.txt"))

##########################################################################################
### Quantify with recountNNLS
##########################################################################################
samp_size=100
rl=150
samps = paste0("sample_", c(paste0("0", 1:9), 10:samp_size))
outdir=paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/ci_tuning/r2/reads/")
bw = paste0(outdir, '/rail/rail-rna_out/coverage_bigwigs/', samps, '.bw')
jx_file = paste0(outdir, '/rail/rail-rna_out/cross_sample_results/junctions.tsv.gz')
table = data.frame(project=rl, run=samps, bigwig_path=bw, rls=rl, paired_end=FALSE)
pheno = processPheno(table)
# load(paste0("/dcl01/leek/data/ta_poc/GencodeV25/list_1p_", rl, ".rda"))
power=0
cores = 20


load("~/downloads/CI_cal_countmat.rda")
load("~/downloads/CI_sim.rda")
count_mat = count_mat/2
rse_sub = rse_tx[match(rownames(count_mat), rownames(rse_tx)),]
bs = assays(rse_sub)$counts
se = assays(rse_sub)$se
score = assays(rse_sub)$scores1
cil = bs-1.96*se
ciu = bs+1.96*se
hit = (cil<=count_mat)*(ciu>=count_mat)

hits = apply(hit, 1, mean)
scatter.smooth(hits~score[,1], col=rgb(1, 0, 0, 0.1), pch=19, ylim=c(0.9, 1))

test = cbind(bs[,1], count_mat[,1], se[,1])






