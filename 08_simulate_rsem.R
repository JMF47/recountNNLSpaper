##### Contains the code to simulate and run analysis based on RSEM estimates of ERR188410

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
		TxL = TxL[TxL$tx_len>150,]
	TxSeq = extractTranscriptSeqs(hg38, TxDb)
		TxSeq = TxSeq[width(TxSeq)>150]
		names(TxSeq) = TxL$tx_name

##########################################################################################
### Estimate the simulation truth from RSEM results on ERR188410
##########################################################################################

rsem = read.table('/dcl01/leek/data/ta_poc/geuvadis/rsem/ERR188410.isoforms.results', header=TRUE)
counts = round(rsem$expected_count)
	mat = match(names(TxSeq), rsem$transcript_id, nomatch=0)
	sim_counts = rep(0, length(TxSeq))
	sim_counts[mat!=0] = counts[mat]
TxL$counts = sim_counts
save(TxL, file='/dcl01/leek/data/ta_poc/geuvadis/sim_rsem.rda')

##########################################################################################
### Simulate the reads using polyester
##########################################################################################
rl=75
paired=TRUE
outdir=paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/rsem_based/", rl, "_", paired*1+1)
dir.create(outdir); setwd(outdir)
writeXStringSet(TxSeq, "tx.fasta")
count_mat = matrix(TxL$counts, ncol=1)
set.seed(rl+paired)
rownames(count_mat) = names(TxSeq)
save(count_mat, file='truth.rda')
polyester::simulate_experiment_countmat('tx.fasta', readmat=count_mat, outdir="reads", paired=paired, readlen=rl)

##########################################################################################
### Quantify
##########################################################################################
rl=75
paired=TRUE
setwd('/dcl01/leek/data/ta_poc/geuvadis/')
samples = "sample_01"
condition = paste0(rl, "_", paired*1+1)
outdir = paste0('/dcl01/leek/data/ta_poc/geuvadis/simulation/rsem_based/', condition)

##########################################################################################
### Rail-RNA
files = paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/rsem_based/", condition, "/reads/", samples, ".fasta")
man = cbind(files, 0, samples)
rail_dir = paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/rsem_based/", condition, "/rail")
dir.create(rail_dir); setwd(rail_dir)
write.table(man, sep="\t", quote=F, col.names=F, row.names=F, file=paste0(rail_dir, "/manifest.txt"))	
files = paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/rsem_based/", condition, "/reads/", samples)
man = cbind(paste0(files, "_1.fasta"), 0, paste0(files, "_2.fasta"), 0, samples)
rail_dir = paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/rsem_based/", condition, "/rail")
dir.create(rail_dir); setwd(rail_dir)
write.table(man, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, file=paste0(rail_dir, "/manifest.txt"))	
setwd(rail_dir)
index="/dcl01/leek/data/HG38_rail/Homo_sapiens/UCSC/hg38/Sequence/"
system2("rail-rna", paste0("go local --deliverables tsv bw -x ", index, "/BowtieIndex/genome ", index, "/Bowtie2Index/genome -m manifest.txt"))

library(recountNNLS); library(recountNNLSdata)
cores=20
samps = "sample_01"
bw = paste0('/dcl01/leek/data/ta_poc/geuvadis/simulation/rsem_based/75_2/rail/rail-rna_out/coverage_bigwigs/', samps, '.bw')
jx_file = paste0('/dcl01/leek/data/ta_poc/geuvadis/simulation/rsem_based/75_2/rail/rail-rna_out/cross_sample_results/junctions.tsv.gz')
table = data.frame(project=rl, run=samps, bigwig_path=bw, rls=rl, paired_end=TRUE)
pheno = processPheno(table)
rse_tx = recountNNLS(pheno, jx_file=jx_file, cores=20)
save(rse_tx, file="/dcl01/leek/data/ta_poc/geuvadis/simulation/rsem_based/75_2/rse_tx_0.rda")

##########################################################################################
### HISAT2-cufflinks (2.0.5-2.2.1)
rl=75
paired=TRUE
outdir=paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/rsem_based/", rl, "_", paired*1+1)
setwd(outdir)
dir.create("hisat2")
dir.create("hisat2/aligned")
dir.create("hisat2/sorted")

s="sample_01"
print(s); flush.console()
outPath = paste0(outdir, "/hisat2/sorted/", s, ".cl.bam")
outBam = paste0(outdir, "/hisat2/sorted/", s, ".cl")
outSam = paste0(outdir, "/hisat2/aligned/", s, ".cl.sam")
intermediateBam = paste0(outdir, "/hisat2/aligned/", s, ".cl.bam")
readPath = paste0(outdir, "/reads/", s)
system2("hisat2", paste0("--dta-cufflinks -f -x /dcl01/leek/data/ta_poc/hisat2-indices/hisat2-hg38-index -1 ", readPath, "_1.fasta", " -2 ", readPath, "_2.fasta ", "-S ", outSam))
system2("samtools", paste0("view -b ", outSam, " > ", intermediateBam))
system2("samtools", paste0("sort ", intermediateBam, " ", outBam))

cufflinks_dir = paste0(outdir, "/hisat2/cufflinks")
system2("mkdir", cufflinks_dir)
cufflinks_dir_default = paste0(cufflinks_dir, "/default")
system2("mkdir", cufflinks_dir_default)
annotation="/dcl01/leek/data/ta_poc/GencodeV25/gencodeV25.coding.gtf"

workingDir = paste0(cufflinks_dir, "/", s)
system2("mkdir", workingDir)
system2("cufflinks", paste0("-m 250 -s 25 --total-hits-norm --no-effective-length-correction --no-length-correction -o ",
	workingDir, " -G ", annotation, " ", outdir, "/hisat2/sorted/", s, ".cl.bam"))

##########################################################################################
### Kallisto (0.43.0)
# PATH=$PATH:/users/jmfu/kallisto_linux-v0.43.0/
rl=75
paired=TRUE
outdir=paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/rsem_based/", rl, "_", paired*1+1)
kallisto_dir = paste0(outdir, "/kallisto")
system2("mkdir", kallisto_dir)
index="/dcl01/leek/data/ta_poc/geuvadis/kallisto/kallisto_index"
s = "sample_01"
system2("kallisto", paste0("quant -i ", index, " -o ", outdir, "/kallisto" , " ", outdir, "/reads/", s, "_1.fasta ", outdir, "/reads/", s, "_2.fasta"))
	
##########################################################################################
### Salmon (0.8.2)
rl=75
paired=TRUE
outdir=paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/rsem_based/", rl, "_", paired*1+1)
salmon_dir = paste0(outdir, "/salmon")
dir.create(salmon_dir)
index="/dcl01/leek/data/ta_poc/GencodeV25/hg38_coding_salmon_ind "
s = "sample_01"

system2("/users/jmfu/Salmon-0.8.2_linux_x86_64/bin/salmon", paste0("quant -i ", index, " -l A -1 ", outdir, "/reads/", s, "_1.fasta -2 ", outdir, "/reads/", s, "_2.fasta -o ", salmon_dir, "/", s))	

#################################################################################
### Compile information
#################################################################################
rm(list=ls())
load('/dcl01/leek/data/ta_poc/geuvadis/simulation/rsem_based/75_2/truth.rda')
load("/dcl01/leek/data/ta_poc/geuvadis/simulation/rsem_based/75_2/rse_tx_0.rda")
ind = match(rownames(count_mat), rownames(rse_tx))
rse_sub = rse_tx[ind,]
tx_info = rowData(rse_sub)
se = assays(rse_sub)$se
score = assays(rse_sub)$score
df = assays(rse_sub)$df

setwd('/dcl01/leek/data/ta_poc/geuvadis/simulation/rsem_based/75_2/')
cl = rtracklayer::import('hisat2/cufflinks/sample_01/transcripts.gtf')
kallisto = read.table('kallisto/abundance.tsv', header=T)
salmon = read.table('salmon/sample_01/quant.sf', header=T)

tx_list = rownames(rse_sub)
out = assays(rse_sub)$counts
kallisto_mat = match(tx_list, kallisto$target_id)
	out = cbind(out, kl = kallisto$est_counts[kallisto_mat])
cl = cl[cl$type=="transcript"]
	cl_mat = match(tx_list, cl$transcript_id)
	cl_cov=as.numeric(cl$cov[cl_mat])
	cl_cts = cl_cov*tx_info$tx_len/75/2
	out = cbind(out, cl_cts)
tx_list = stringr::str_replace(tx_list, "_PAR_Y", "")
salmon_mat = match(tx_list, salmon$Name)
	out = cbind(out, sl = salmon$NumReads[salmon_mat])
	out = data.frame(out)
save(out, count_mat, se, score, df, file="~/rsem_based_0.rda")

