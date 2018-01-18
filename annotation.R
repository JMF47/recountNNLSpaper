rm(list=ls())
library(rtracklayer); library(exomeCopy)
library(Biostrings); library(exomeCopy); library(GenomicFeatures); library(BSgenome.Hsapiens.UCSC.hg38)

### Download GencodeV25 annotation
download.file("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz",
	destfile="/dcl01/leek/data/ta_poc/GencodeV25/gencodeV25.gff3.gz")
system2("gzip", "-d /dcl01/leek/data/ta_poc/GencodeV25/gencodeV25.gff3.gz")

### Extract protein_coding genes
gff_all = import.gff("/dcl01/leek/data/ta_poc/GencodeV25/gencodeV25.gff3")
gff_coding = gff_all[which(gff_all$gene_type=="protein_coding")]
export(gff_coding, con="/dcl01/leek/data/ta_poc/GencodeV25/gencodeV25.coding.gff3")

###### Calculating features for recountNNLS
### Disjoin exons
gff_coding = import.gff("/dcl01/leek/data/ta_poc/GencodeV25/gencodeV25.coding.gff3")
gff_exons = gff_coding[gff_coding$type=="exon"]
gff_dj = disjoin(gff_exons, ignore.strand=TRUE)
export(gff_dj, con="/dcl01/leek/data/ta_poc/GencodeV25/bins_dj.bed", format="bed")

### Further subdivide the exonic segments
### So that no segment is more than twice the length of the read length of experiment
sGR = function(i, subset, subsize){
	print(i); flush.console()
	return(subdivideGRanges(subset[i], subsize=subsize))
}

gff_dj = import.bed("/dcl01/leek/data/ta_poc/GencodeV25/bins_dj.bed")
rls = c(37, 50, 75, 100, 150)
for(rl in rls){
	gff_nodivide = gff_dj[width(gff_dj)<(rl*2)]
	values(gff_nodivide) = NULL
	gff_divide = gff_dj[width(gff_dj)>=(rl*2)]
	subset_divided = mclapply(1:length(gff_divide), sGR, gff_divide, rl, mc.cores=detectCores())
	subset_grl = GRangesList(unlist(subset_divided))
	subset_ul = unlist(subset_grl)
	total = c(gff_nodivide, subset_ul)
	total_clean = sortSeqlevels(total)
	total_clean = sort(total_clean, ignore.strand=TRUE)
	export.bed(total_clean, con=paste0("/dcl01/leek/data/ta_poc/GencodeV25/bins_", rl, ".bed"), format="bed")
}

### Extract the introns
getIntrons = function(chr, gff){
	print(chr); flush.console()
	gff_sub = gff[seqnames(gff)==chr]
	introns_start_plus = end(gff_sub)[-length(gff_sub)]+1
	introns_end_plus = start(gff_sub)[-1]-1
	introns_end_minus = start(gff_sub)[-length(gff_sub)]-1
	introns_start_minus = end(gff_sub)[-1]+1
	intron_info = cbind(c(1, introns_start_plus), c(1, introns_end_plus), c(1, introns_start_minus), c(1, introns_end_minus))
	intron_use = intron_info[,1:2]; intron_use[which(strand(gff_sub)=="-"),] = intron_info[which(strand(gff_sub)=="-"),3:4]
	txs = unique(gff_sub$transcript_id)
	rm_ind = match(txs, gff_sub$transcript_id)
	output = GRanges(seqnames=chr, IRanges(intron_use[-rm_ind,1], intron_use[-rm_ind,2]), strand=strand(gff_sub)[-rm_ind])
	return(output)
}
gff_exons = gff_coding[gff_coding$type=="exon"]
gff_exons_sub = gff_exons[seqnames(gff_exons)=="chr1"]

gff_jx = unique(unlist(GRangesList(sapply(paste0("chr", c(1:22, "X", "Y")), getIntrons, gff_exons))))
save(gff_jx, file="/dcl01/leek/data/ta_poc/GencodeV25/gff_jx.rda")
