module load python/2.6.6
#### Using Rail to get P
library(GenomicRanges); library(rtracklayer); library(stringr); library(Biostrings); library(Rsamtools); 
library(Biostrings); library(GenomicFeatures); library(BSgenome.Hsapiens.UCSC.hg38)
setwd("/dcl01/leek/data/ta_poc/")
### Get all the transcript sequences
library(BSgenome.Hsapiens.UCSC.hg38)
hg38 = BSgenome.Hsapiens.UCSC.hg38
TxDb = makeTxDbFromGFF("GencodeV25/gencodeV25.coding.gff3")
	TxL = transcriptLengths(TxDb)
	TxSeq = extractTranscriptSeqs(hg38, TxDb)
	names(TxSeq) = TxL$tx_name

#### Start with
# ENSG00000229571.6
# ENSG00000277058.2
# ENSG00000204505.4

gene_id = 'ENSG00000188257.10'
setwd('/dcl01/leek/data/ta_poc/GencodeV25/pileup_150/rail')
TxSeq = import(paste0(gene_id, ".fasta"))

writeReads = function(sequence, readLen, directory){
	suppressWarnings(dir.create(directory))
	indz = 1:(width(as.character(sequence))-readLen+1)
	subReads = DNAStringSet(str_sub(sequence, indz, indz+readLen-1))
	fileOut = paste0(names(sequence), ".fa")
	names(subReads) = paste0(names(sequence), "-", indz)
	writeXStringSet(subReads, file=paste0(directory, "/", fileOut), format="fasta")
	return(fileOut)
}
.processSample = function(sampleFile, grl, bins, verbose=TRUE){
      if(verbose==TRUE)
            message("Processing sample ", sampleFile)

      cov_rle = rtracklayer::import(sampleFile, as = 'RleList')
      ## To ensure consistency when some samples have chrs dropped from no mapped reads
      cov_rle_matched = cov_rle[match(names(grl), names(cov_rle), nomatch=0)]
      grl_keep = grl[match(names(cov_rle_matched), names(grl), nomatch=0)]
      cov_binned = sapply(names(grl_keep), .processChr, cov_rle_matched, grl_keep)
            if(class(cov_binned)=="list")
            	cov_binned = do.call(c, cov_binned)
      id = queryHits(GenomicRanges::findOverlaps(bins, unlist(grl_keep), type="equal"))

      cov_out = rep(0, length(bins))
      cov_out[id[!is.na(id)]] = cov_binned[!is.na(id)]
      return(cov_out)
}
.processChr = function(chr, rle_cov, grl_bins){
      sum(Views(rle_cov[[chr]], rtracklayer::ranges(grl_bins[[chr]])))
}

####################################################################################################
##### Ex ############################################################
txs = names(TxSeq)
tx_seqs = TxSeq[txs]

readLens = 150
directory = paste0("/dcl01/leek/data/ta_poc/GencodeV25/pileup_150/rail/")
filesOut = NULL
for(i in 1:length(tx_seqs)){
	filesOut = c(filesOut, writeReads(tx_seqs[i], readLens, paste0(directory, gene_id)))
}

files = paste0(directory, gene_id, "/", unlist(filesOut))
man = cbind(files, 0, txs)
rail_dir = paste0(directory, gene_id)
setwd(rail_dir)
write.table(man, sep="\t", quote=F, col.names=F, row.names=F, file="manifest.txt")
index="/dcl01/leek/data/HG38_rail/Homo_sapiens/UCSC/hg38/Sequence/"
system2("rail-rna", paste0("go local --deliverables tsv bw -x ", index, "/BowtieIndex/genome ", index, "/Bowtie2Index/genome -m ", rail_dir, "/manifest.txt"))

bw_files = list.files(paste0(rail_dir, "/rail-rna_out/coverage_bigwigs"))
	bw_files = bw_files[match(paste0(txs, ".bw"), bw_files)]

data(bins_150, package="recountNNLSdata")
gff_ex = bins_150

### Count to features (exons)
grl = GenomicRanges::split(gff_ex, GenomicRanges::seqnames(gff_ex))
list_totCov= lapply(paste0(rail_dir, "/rail-rna_out/coverage_bigwigs/", bw_files), .processSample, grl, gff_ex)
F_ex = do.call(cbind, list_totCov)
	rownames(F_ex) = paste0('e', 1:dim(F_ex)[1])
	F_ex = F_ex[rowSums(F_ex)>0,]
	F_ex = t(t(F_ex)/(width(tx_seqs)-150+1))
	colnames(F_ex) = txs


##### Jx ############################################################
jx_files = list.files(paste0(rail_dir, "/rail-rna_out/cross_sample_results"))
jx_file = paste0(rail_dir, "/rail-rna_out/cross_sample_results/junctions.tsv.gz")
F_jx = getJxCounts(jx_file)
	F_jx = t(t(F_jx)/(width(tx_seqs)-150+1))
F = rbind(F_ex/150, F_jx[,match(colnames(F_ex), colnames(F_jx))])
colnames(F) = txs


##########################################################################################
base_dir = "/dcl01/leek/data/ta_poc/geuvadis"
setwd(base_dir)
rl = 37
paired = 1
condition = paste0(rl, "_", paired)
samples = paste0("sample_0", 1)
outdir = paste0("simulation/", condition)
dir.create(outdir)
setwd(outdir)

# samps = c("sample_01", "sample_02")
samps = "sample_01"
bw = paste0('/dcl01/leek/data/ta_poc/geuvadis/simulation/', condition, '/rail/rail-rna_out/coverage_bigwigs/', samps, '.unique.bw')
jx_file = paste0('/dcl01/leek/data/ta_poc/geuvadis/simulation/', condition, '/rail/rail-rna_out/cross_sample_results/junctions.tsv.gz')
table = data.frame(project=rl, run=samps, bigwig_path=bw, rls=rl, paired_end=(paired==2))
pheno = processPheno(table)
counts = getCounts(pheno, jx_file)

counts_sub = counts[match(rownames(F_ex), rownames(counts)),]
F_ex37 = F_ex/37
test = nnls::nnls(F_ex37, counts_sub)


test = sapply(counts_sub, .lnnls, F_ex/37)


.lnnls = function(data, matrix, boot=100){
      mat = matrix(apply(matrix, 2, as.numeric), nrow=dim(matrix)[1])
            # penalty = 1+1*(data==0)
            # mat = mat*penalty
            # data = data*penalty
      mod=nnls::nnls(mat, data)
      return(mod)
}





