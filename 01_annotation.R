##### Contains the code to download and process the necessary
##### annotation and reference fasta

##########################################################################################
### Obtain the annotation and reference fasta
##########################################################################################
rm(list=ls())
library(polyester); library(rtracklayer); library(Biostrings); library(stringr)

system2("wget", "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz")
system2("gunzip", "hg38.fa.gz")
system2("gffread", "-w sim.fa -g hg38.fa gencodeV25.sim.gtf")
hg38 = import("hg38.fa", format="fasta")
hg38 = hg38[names(hg38) %in% paste0("chr", c(1:22, "X", "Y", "M"))]
	export(hg38, con="hg38.fa", format="fasta")
hg38_1_14 = hg38[names(hg38) %in% c("chr1", "chr14")]
	export(hg38_1_14, con="hg38_1_14.fa", format="fasta")

### Feel free to change to a different directory
setwd("~")
download.file("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz",
	destfile="gencodeV25.gff3.gz")
system2("gzip", "-d gencodeV25.gff3.gz")

### Extract protein_coding genes
gff_all = import.gff("gencodeV25.gff3")
gff_coding = gff_all[which(gff_all$gene_type=="protein_coding")]
	export(gff_coding, con="gencodeV25.coding.gtf")

##########################################################################################
### Subset to coding transcripts on chr1 and chr14
##########################################################################################
txs_sim = gff_coding[seqnames(gff_coding) %in% c("chr1", "chr14")]
	export(txs_sim, con='gencodeV25.sim.gtf', format='gtf')