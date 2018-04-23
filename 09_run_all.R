rm(list=ls())
library(GenomicRanges); library(rtracklayer); library(stringr); library(Biostrings); library(Rsamtools)
library(SummarizedExperiment); library(recountNNLSdata); library(recountNNLS)

url_table <- recount::recount_url
	unique_ids = unique(url_table$project)
unique_ids = as.character(unique_ids[! unique_ids %in% c("TCGA", "SRP012682")])

for(unique_id in unique_ids){
    message(which(unique_ids==unique_id))
    out_dir = paste0("/dcl01/leek/data/ta_poc/recount_out/rse_new/", unique_id)
    if(dir.exists(out_dir)==F){
    	dir.create(out_dir)
		setwd(out_dir)
		pheno = processPheno(unique_id, local=TRUE)
		rse_tx = recountNNLS(pheno, cores=20)	
		save(rse_tx, file="rse_tx.RData")	
		rm(rse_tx, pheno, out_dir)
    }
}

### Spot-check
rm(list=ls())
library(GenomicRanges); library(rtracklayer); library(stringr); library(Biostrings); library(Rsamtools)
library(SummarizedExperiment); library(recountNNLSdata); library(recountNNLS)

url_table <- recount::recount_url
	unique_ids = unique(url_table$project)
unique_ids = as.character(unique_ids[! unique_ids %in% c("TCGA", "SRP012682")])

out = NULL
for(unique_id in unique_ids[736:length(unique_ids)]){
	message(which(unique_ids == unique_id))
	out_dir = paste0("/dcl01/leek/data/ta_poc/recount_out/rse_new/", unique_id)
	out_file = paste0(out_dir, "/rse_tx.RData")
	pheno = processPheno(unique_id, local=TRUE)
	load(out_file)
	pheno_count = pheno$mapped_read_count/(pheno$paired_end+1)
	nnls_count = apply(assays(rse_tx)$fragments, 2, sum, na.rm=TRUE)
	out = rbind(out, cbind(pheno_count, nnls_count))
}

