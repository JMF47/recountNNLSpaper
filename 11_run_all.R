#### Pipeline to estimate all recount samples
### (minus TCGA and GTEX)

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
