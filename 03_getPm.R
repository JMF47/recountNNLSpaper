##### Contains the code to create probability matrices
rm(list=ls())
library(GenomicRanges); library(rtracklayer); library(stringr); library(Biostrings); library(Rsamtools); 
library(Biostrings); library(GenomicFeatures); library(BSgenome.Hsapiens.UCSC.hg38)

##################################################################################
### Create the bins for each read length category
##################################################################################
library(rtracklayer)
gff = import.gff("/dcl01/leek/data/ta_poc/GencodeV25/gencodeV25.coding.gtf")
gff_exons = gff[gff$type=="exon"]
gff_jx = unique(unlist(GRangesList(mclapply(paste0("chr", c(1:22, "X", "Y")), getIntrons, gff_exons, mc.cores=20))))
	save(gff_jx, file="/dcl01/leek/data/ta_poc/GencodeV25/gff_jx.rda")
gff_dj = disjoin(gff_exons, ignore.strand=TRUE)
export(gff_dj, con="/dcl01/leek/data/ta_poc/GencodeV25/bins_dj.bed", format="bed")

gff_ex = import.bed("/dcl01/leek/data/ta_poc/GencodeV25/bins_dj.bed")

sGR = function(i, subset, subsize){
	print(i); flush.console()
	return(exomeCopy::subdivideGRanges(subset[i], subsize=subsize))
}

rls = c(37, 50, 75, 100, 150)
for(rl in rls){
	gff_nodivide = gff_ex[width(gff_ex)<(rl*2)]
		values(gff_nodivide) = NULL
	gff_divide = gff_ex[width(gff_ex)>=(rl*2)]
	subset_divided = mclapply(1:length(gff_divide), sGR, gff_divide, rl, mc.cores=20)
	subset_grl = GRangesList(unlist(subset_divided))
	subset_ul = unlist(subset_grl)
	total = c(gff_nodivide, subset_ul)
	total_clean = sortSeqlevels(total)
	total_clean = sort(total_clean, ignore.strand=T)
	export.bed(total_clean, con=paste0("/dcl01/leek/data/ta_poc/GencodeV25/bins_", rl, ".bed"), format="bed")
}

##################################################################################
### Functions for matrix creation
##################################################################################
getSequence = function(transcript_id, TxSeq){
	return(TxSeq[names(TxSeq)==transcript_id])
}
writeReads = function(sequence, readLen, directory){
	suppressWarnings(dir.create(directory))
	indz = 1:max(1, (width(sequence)-readLen+1))
	subReads = DNAStringSet(str_sub(sequence, indz, indz+readLen-1))
	fileOut = paste0(names(sequence), ".fa")
	names(subReads) = paste0(names(sequence), "-", indz)
	writeXStringSet(subReads, file=paste0(directory, "/", fileOut), format="fasta")
	return(fileOut)
}
alignReads = function(filename, aligner="hisat2", directory, gene, ref){
	indexPath = "/dcl01/leek/data/ta_poc/hisat2-indices/hisat2-hg38-index"
	pathIn = paste0(directory, gene,"/", filename)
	pathOutSam = paste0(pathIn, ".sam")
	pathOutBam = paste0(pathOutSam, ".bam")
	pathOutSplit = paste0(pathOutBam, ".split")
	system2("hisat2", paste0("-x ", indexPath, " -f -U ", pathIn," -S ", pathOutSam))
	system2("samtools", paste0("sort ", pathOutSam, " ", pathOutSam))
	system2("samtools", paste0("index ", pathOutBam))
	system2("regtools", paste("junctions extract -i 10 -o", pathOutSplit, pathOutBam))
	return(pathOutBam)	
}
getBamGrJunctions = function(bamPath){
	splitPath = paste0(bamPath, ".split")
	if(file.size(splitPath)>0){
		junction_info = read.table(splitPath, colClasses = c("character", "integer", "integer", "NULL", "integer", "NULL", "NULL", "NULL", "NULL", "NULL", "character", "NULL"))
		if(str_detect(as.character(junction_info[1,1]), "chr")==FALSE)
			junction_info[,1] = paste0("chr", str_replace(as.character(junction_info[,1]), "T", ""))
		blocks = str_split_fixed(junction_info[,5], ",", n=2)
		junction_gr = GRanges(seqnames=junction_info[,1], IRanges(junction_info[,2]+1+as.numeric(blocks[,1]), junction_info[,3]-as.numeric(blocks[,2])), score = junction_info[,4])
		return(junction_gr)
	}
	return(NULL)
}
evalPileupJunctions = function(gr, gff_jx){
	ol_jx = rep(0, length(gff_jx))
	if(length(gr)>0){
		ol = findOverlaps(gr, gff_jx, type="equal")
		ol_jx[subjectHits(ol)] = gr$score[queryHits(ol)]
		return(ol_jx)
	}
	return(ol_jx)
}
evalPileupExons = function(bamPath, gff_ex){
	ol_ex = rep(0, length(gff_ex))
	system2("samtools", paste0("depth ", bamPath, " > ", bamPath, ".etab"))
	count_file = paste0(bamPath, ".etab")
	if(file.size(count_file)>0){
		covinfo = read.table(count_file)
		covinfo = covinfo[covinfo[,3]!=0,]
		if(str_detect(as.character(covinfo[1,1]), "chr")==FALSE)
			covinfo[,1] = paste0("chr", str_replace(as.character(covinfo[,1]), "T", ""))
		covinfo_gr = GRanges(covinfo[,1], IRanges(covinfo[,2], covinfo[,2]), depth = covinfo[,3])
		ol = findOverlaps(covinfo_gr, gff_ex)
		if(length(ol)>0){
			depth_info = by(covinfo_gr$depth[queryHits(ol)], subjectHits(ol), sum)
			ol_ex[as.numeric(names(depth_info))] = as.numeric(depth_info)	
		}	
	}
	return(ol_ex)
}
writeF = function(gene_id, directory, readLens, gff_ex, gff_jx){
	print(which(genes==gene_id)); flush.console()
	tabOut = paste0(directory, "Fs/", gene_id, ".tab")
	if(!file.exists(tabOut)){
		file.create(tabOut)
		txs = TxL$tx_name[TxL$gene_id==gene_id]

		# Obtain sequence
		sequences_txs = sapply(txs, getSequence, TxSeq)
		if(length(sequences_txs)>0){
			# Write sequence
			filesOut = lapply(sequences_txs, writeReads, readLens, paste0(directory, gene_id))
			# Align sequence
			bamPaths = lapply(filesOut, alignReads, directory=directory, gene=gene_id)
			# GRs
			grs_junctions = lapply(bamPaths, getBamGrJunctions)
			counts_junctions = do.call(cbind, lapply(grs_junctions, evalPileupJunctions, gff_jx))
			counts_exons = do.call(cbind, lapply(bamPaths, evalPileupExons, gff_ex))
				# Check if any reads mapped at all
				read_counts = apply(counts_exons, 2, sum) + apply(counts_junctions, 2, sum)
				keep_ind = which(read_counts>0)
			# If atleast 1 transcript made it through the alignment + counting
			if(length(keep_ind)>0){
				rownames(counts_exons) = paste0("e", 1:length(gff_ex))
				rownames(counts_junctions) = paste0("i", 1:length(gff_jx))
				C = rbind(counts_exons[rowSums(counts_exons)>0,, drop=FALSE]/readLens, 
					counts_junctions[rowSums(counts_junctions)>0,,drop=FALSE])
				C = C[,keep_ind, drop=FALSE]
				if(length(C)>0){
					colnames(C) = txs[keep_ind]
					wids = unlist(lapply(sequences_txs, function(x) max(readLens, nchar(x))))
					F = t(t(C)/(unlist(wids[keep_ind])-readLens+1))
					write.table(F, quote=FALSE, col.names=T, row.names=T, file=tabOut)
				}
			}
			system2("rm", paste0("-r ", directory, gene_id))
		}
	}
	return(0)
}

##################################################################################
### Obtain matrices for each read length
##################################################################################
setwd("/dcl01/leek/data/ta_poc/")
### Get all the transcript sequences
library(BSgenome.Hsapiens.UCSC.hg38)
hg38 = BSgenome.Hsapiens.UCSC.hg38
TxDb = makeTxDbFromGFF("GencodeV25/gencodeV25.coding.gff3")
	TxL = transcriptLengths(TxDb)
	TxSeq = extractTranscriptSeqs(hg38, TxDb)
	names(TxSeq) = TxL$tx_name
load("GencodeV25/gff_jx.rda")
genes = unique(TxL$gene_id)
rls = c(37, 50, 75, 100, 150)

load('/dcl01/leek/data/ta_poc/GencodeV25/gff_jx.rda')
for(readLens in rls){
	directory = paste0("GencodeV25/pileup_", readLens, "/")
	gff_ex = import.bed(paste0("GencodeV25/bins_", readLens, ".bed"))
	diagnostic = mclapply(genes, writeF, directory, readLens, gff_ex, gff_jx, mc.cores=20)	
}

##########################################################################################
### Compile matrices into rdas
##########################################################################################
library(stringr); library(parallel)
getEmissionHelper = function(file, rl){
	message(which(files==file))
	emission_matrix = NULL
	if(file.size(file)>0){
		emission_matrix = read.table(file)
	}
	return(emission_matrix)
}

rls = c(37, 50, 75, 100, 150)
for(rl in rls){
	setwd(paste0("/dcl01/leek/data/ta_poc/GencodeV25/pileup_", rl, "/Fs"))
	files = list.files()
	matrix_list = mclapply(files, getEmissionHelper, rl, mc.cores=20)
	names(matrix_list) = stringr::str_replace(files, ".tab", "")
	assign(paste0("matrix_", rl), matrix_list)
	save(list = paste0('matrix'_, rl), 
		file=paste0('/dcl01/leek/data/ta_poc/GencodeV25/matrix_', rl, '.rda'))
}




