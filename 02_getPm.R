##### Contains the code create probability matrices
rm(list=ls())
library(GenomicRanges); library(rtracklayer); library(stringr); library(Biostrings); library(Rsamtools); 
library(Biostrings); library(exomeCopy); library(GenomicFeatures); library(BSgenome.Hsapiens.UCSC.hg38)

########################################################################
### functions for matrix creation
getSequence = function(transcript_id, TxSeq){
	return(TxSeq[names(TxSeq)==transcript_id])
}
writeReads = function(sequence, readLen, directory){
	suppressWarnings(dir.create(directory))
	indz = 1:(width(sequence)-readLen+1)
	subReads = DNAStringSet(str_sub(sequence, indz, indz+readLen-1))
	fileOut = paste0(names(sequence), ".fa")
	names(subReads) = paste0(names(sequence), "-", indz)
	writeXStringSet(subReads, file=paste0(directory, "/", fileOut), format="fasta")
	return(fileOut)
}
alignReads = function(filename, aligner="hisat2", directory, gene, ref){
	indexPath = "hisat2-index-hg38/hisat2-ind"
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
		
		minLen = sapply(sequences_txs, width)>=readLens
			rm_ind = which(minLen==FALSE)
			if(length(rm_ind)>0){
				sequences_txs = sequences_txs[-rm_ind]
				txs = txs[-rm_ind]		
			}
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
			if(length(keep_ind)>0){
				rownames(counts_exons) = paste0("e", 1:length(gff_ex))
				rownames(counts_junctions) = paste0("i", 1:length(gff_jx))
				C = rbind(counts_exons[rowSums(counts_exons)>0,, drop=FALSE], counts_junctions[rowSums(counts_junctions)>0,,drop=FALSE])
				C = C[,keep_ind, drop=FALSE]
				if(length(C)>0){
					colnames(C) = txs[keep_ind]
					wids = unlist(lapply(sequences_txs, nchar))
					F = t(t(C)/(unlist(wids[keep_ind])-readLens+1))
					write.table(F, quote=FALSE, col.names=T, row.names=T, file=tabOut)
				}
			}
			system2("rm", paste0("-r ", directory, gene_id))
		}
	}
	return(0)
}
########################################################################

### Get all the transcript sequences
hg38 = BSgenome.Hsapiens.UCSC.hg38
TxDb = makeTxDbFromGFF("GencodeV25/gencodeV25.coding.gff3")
	TxL = transcriptLengths(TxDb)
	TxSeq = extractTranscriptSeqs(hg38, TxDb)
	names(TxSeq) = TxL$tx_name
load("GencodeV25/gff_jx.rda")
genes = unique(TxL$gene_id)
rls = c(37, 50, 75, 100, 150)

for(readLens in rls){
	directory = paste0("GencodeV25/pileup_", readLens, "/")
	gff_ex = import.bed(paste0("GencodeV25/bins_", readLens, ".bed"))
	diagnostic = mclapply(genes, writeF, directory, readLens, gff_ex, gff_jx, mc.cores=10)	
}

directory = "GencodeV25/pileup_50/"
	gff_ex = import.bed(paste0("GencodeV25/bins_", readLens, ".bed")
	diagnostic = mclapply(genes, writeF, directory, readLens, gff_ex, gff_jx, mc.cores=20)	

##########################################################################################
### create rda
##########################################################################################
library(stringr); library(parallel)
getEmissionHelper = function(file, rl){
	message(which(files==file))
	emission_matrix = NULL
	if(file.size(file)>0){
		emission_matrix = read.table(file)
		hasE = str_detect(rownames(emission_matrix), "e")
		emission_matrix[hasE,] = emission_matrix[hasE,]/rl
	}
	return(emission_matrix)
}

rl = 75
setwd(paste0("GencodeV25/pileup_", rl, "/Fs"))
files = list.files()
matrix_list = mclapply(files, getEmissionHelper, rl, mc.cores=10)
names(matrix_list) = str_replace(files, ".tab", "")
matrix_75 = matrix_list

### Rescale the exonic feature probabilities to be on read scale
for(i in 1:length(files)){
	print(i); flush.console()
	f = files[i]
	emission_matrix = NULL
	if(file.size(f)>0){
		emission_matrix = read.table(f)
		hasE = str_detect(rownames(emission_matrix), "e")
		emission_matrix[hasE,] = emission_matrix[hasE,]/rl
		matrix_list[[i]] = emission_matrix
	}
}
names(matrix_list) = stringr::str_replace(files, ".tab", "")
save(matrix_list, file='GencodeV25/matrix_175.rda')