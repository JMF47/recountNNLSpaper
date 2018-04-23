#### Pipeline for GTEX samples

##########################################################################################
## Break apart junction file into smaller pieces
## The bed file is in the same order as the junction coverage file
##########################################################################################

##########################################################################################
## Read in where the junctions are
##########################################################################################
rm(list=ls())
library(rtracklayer)
raw=import('/dcl01/leek/data/ta_poc/SRP012682.junction_id_with_transcripts.bed.gz')
end(raw) = end(raw)+1
data(gff_jx, package="recountNNLSdata")
ol = findOverlaps(gff_jx, raw, type="equal")

##########################################################################################
## Read in raw counts
## Filter to only annotated junctions
## Create matrix of only annotated counts
##########################################################################################
tab = read.table('/dcl01/leek/data/ta_poc/SRP012682.junction_coverage.tsv.gz')
tab_sub = tab[subjectHits(ol),]
rm(tab); gc()
samps = read.table("/dcl01/leek/data/ta_poc/sample_ids.tsv")
samps_gtex = samps[samps[,2]=="SRP012682",]
cts = matrix(0, ncol=dim(samps_gtex)[1], nrow=length(ol))

new_out = NULL
inds = NULL
for(i in 1:length(ol)){
	message(i)
	info = tab_sub[i,]
	subjs = unlist(stringr::str_split(info[,2], ","))
	jxcts = as.numeric(unlist(stringr::str_split(info[,3], ",")))
	new = rep(0, dim(meta_use)[1])
		new[match(subjs, colnames(cts))] = jxcts
	new_out = rbind(new_out, new)
	inds = c(inds, i)
	if(i %% 2000 == 0){
		cts[inds, ] = new_out
		new_out = NULL
		inds = NULL
	}
}
cts[inds, ] = new_out
id_map = read.table('https://jhubiostatistics.shinyapps.io/recount/sample_ids.tsv')
colnames(cts) = samps_gtex[,1]
mat = match(colnames(cts), id_map[,1])
colnames(cts) = id_map[mat,3]

##########################################################################################
### Subset counts into chunks of at most 300 samples
##########################################################################################
meta_use <- recountNNLS::processPheno("SRP012682")
ids = paste0(meta_use$project, "-", meta_use$rls_group, "-", 1)
unique_ids = unique(ids)
id_counts = by(rep(1, length(ids)), ids, sum)
for(id in names(id_counts[id_counts>300])){
	ind = which(ids==id)
	replace = rep(1:ceiling(length(ind)/300), each=300)[1:length(ind)]
	ids[ind] = str_replace(ids[ind], "-[0-9]*$", paste0("-", replace))		
}
meta_use$unique_ids = ids

unique_ids = unique(ids)
setwd("/dcl01/leek/data/ta_poc/rse_jx")
for(i in 1:length(unique_ids)){
	message(i)
	unique_id = unique_ids[i]
	meta_sub = meta_use[meta_use$id == unique_id,]
	mat = match(meta_sub$run, colnames(cts))
	cts_sub = cts[,mat]
	rse_jx = SummarizedExperiment::SummarizedExperiment(assays=list(counts=cts_sub), 
		rowRanges=gff_jx[queryHits(ol)], colData=meta_sub)
	file_out = paste0(unique_id, ".rse_jx.rda")
	save(rse_jx, file=file_out)
	rm(rse_jx, meta_sub, mat, cts_sub)
}

##########################################################################################
### Run algorithm on the chunks of samples
##########################################################################################
rm(list=ls())
library(GenomicRanges); library(rtracklayer); library(stringr); library(Biostrings); library(Rsamtools)
library(SummarizedExperiment); library(recountNNLSdata); library(recountNNLS)
setwd("/dcl01/leek/data/ta_poc/rse_jx")
unique_ids = list.files(pattern="SRP012682")
pheno_full = processPheno("SRP012682", local=TRUE)
data(gff_jx, package="recountNNLSdata")

dir.create("/dcl01/leek/data/ta_poc/recount_out/rse_new/SRP012682/")
for(unique_id in unique_ids){
    message(which(unique_ids==unique_id))
    out_dir = paste0("/dcl01/leek/data/ta_poc/recount_out/rse_new/SRP012682/", unique_id)
    if(dir.exists(out_dir)==F){
    	dir.create(out_dir)
		load(unique_id)
		rl = as.numeric(str_replace_all(str_extract(unique_id, "-[0-9]*-"), "-", ''))
		pheno = pheno_full[pheno_full$run %in% colnames(rse_jx),]
		counts_ex = getExCounts(pheno, cores = 1)
      	ol = findOverlaps(gff_jx, SummarizedExperiment::rowRanges(rse_jx), type="equal")
        counts_jx = SummarizedExperiment::assays(rse_jx[, match(pheno$run, colnames(rse_jx))])$counts
		rownames(counts_jx) = paste0("i", queryHits(ol))
		samps = intersect(colnames(counts_ex), colnames(counts_jx))
		counts = rbind(counts_ex[,samps,drop=FALSE], counts_jx[,samps,drop=FALSE])
		rse_tx = processReadLength(rl, pheno, NULL, 20, counts)
		save(rse_tx, file=paste0(out_dir, "/rse_tx.rda"))	
		rm(rse_tx, pheno, out_dir)
    }
}

##########################################################################################
### Merge chunks and then separate by tissue
##########################################################################################
rm(list=ls())
library(GenomicRanges); library(rtracklayer); library(stringr); library(Biostrings); library(Rsamtools)
library(SummarizedExperiment); library(recountNNLSdata); library(recountNNLS)
setwd("/dcl01/leek/data/ta_poc/rse_jx")
unique_ids = list.files(pattern="SRP012682")
out_list = list()
for(i in 1:length(unique_ids)){
	message(i)
	unique_id = unique_ids[i]
	file = paste0("/dcl01/leek/data/ta_poc/recount_out/rse_new/SRP012682/", unique_id, "/rse_tx.rda")
	load(file)
	out_list[[i]] = rse_tx
}
rse_tx = do.call(SummarizedExperiment::cbind, out_list)
save(rse_tx, file='/dcl01/leek/data/ta_poc/recount_out/rse_new/SRP012682/rse_tx.RData')

fix_tissue <- function(rse) {
    ## Fix 5 samples that are missing smts info but have smtsd info
    colData(rse)$smts[colData(rse)$smts == '' & colData(rse)$smtsd == 'Esophagus - Mucosa'] <-  'Esophagus'
    colData(rse)$smts[colData(rse)$smts == '' & colData(rse)$smtsd == 'Skin - Sun Exposed (Lower leg)'] <-  'Skin'
    colData(rse)$smts[colData(rse)$smts == '' & colData(rse)$smtsd == 'Stomach'] <-  'Stomach'
    
    ## Print number of samples by tissue
    print(sort(table(colData(rse)$smts)))
    
    ## Finish
    return(rse)
}

split_by_tissue <- function(rse, type) {
    ## Split data by tissue
    tissues <- unique(colData(rse)$smts)
    
    message(paste(Sys.time(), 'splitting', type, 'level information by tissue'))
    rse_split <- lapply(tissues, function(tissue) {
        subset(rse, select = colData(rse)$smts == tissue)
    })
    names(rse_split) <- tissues
    
    ## Run some checks
    stopifnot(sum(sapply(rse_split, ncol)) == ncol(rse))
    stopifnot(all(sort(sapply(rse_split, ncol)) - sort(table(colData(rse)$smts)) == 0))
    
    ## For file names
    tissues <- gsub(' ', '_', tolower(tissues))
    
    res <- mapply(function(rse_tissue, tissue, type) {
        
        message(paste(Sys.time(), 'saving', type,
            'level information for tissue', tissue))
        
        rse_file <- paste0('/dcl01/leek/data/ta_poc/recount_out/rse_new/SRP012682/rse_', type, '_', tissue, '.Rdata')
        
        ## Make sure the variable name is rse_exon or rse_gene
        varname <- paste0('rse_', type)
        assign(varname, rse_tissue)
        
        ## Save
        save(list = varname, file = rse_file)
        return(rse_file)
    }, rse_split, tissues, MoreArgs = list(type = type), SIMPLIFY = FALSE)
}

message(paste(Sys.time(), 'splitting by tissue at the exon level'))
split_by_tissue(fix_tissue(rse_tx), 'tx')
