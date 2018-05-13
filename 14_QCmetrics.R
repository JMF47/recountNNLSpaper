### Spot-check
rm(list=ls())
library(GenomicRanges); library(rtracklayer); library(stringr); library(Biostrings); library(Rsamtools)
library(SummarizedExperiment); library(recountNNLSdata); library(recountNNLS)

url_table <- recount::recount_url
	unique_ids = unique(url_table$project)
unique_ids = as.character(unique_ids[! unique_ids %in% c("TCGA", "SRP012682")])

# out = NULL
overall_map = NULL
for(unique_id in unique_ids){
	message(which(unique_ids == unique_id))
	out_dir = paste0("/dcl01/leek/data/ta_poc/recount_out/rse_new/", unique_id)
	out_file = paste0(out_dir, "/rse_tx.RData")
	load(out_file)
	pheno = colData(rse_tx)
	total_reads = pheno$reads_downloaded/(pheno$paired_end+1)
	mapped_reads = pheno$mapped_read_count/(pheno$paired_end+1)
	tx_reads = apply(assays(rse_tx)$fragments, 2, sum, na.rm=TRUE)
	time = pheno$biosample_submission_date
	info = cbind(pheno[,c(1:4, 6, 8:9)], time, pheno$rls, pheno$rls_group, total_reads, mapped_reads, tx_reads)
	overall_map = rbind(overall_map, info)
}

map_perc = overall_map[,12]/overall_map[,11]
tx_perc = overall_map[,13]/overall_map[,12]
num_fail_tx = by(tx_perc, overall_map$project, function(x) sum(x<0.5))
num_fail_map = by(map_perc, overall_map$project, function(x) sum(x<0.8))
num_sample = by(rep(1, length(map_perc)), overall_map$project, sum)
cbind(num_fail_tx, num_fail_map, num_sample)

quantile(tx_perc, seq(0, 1, by=0.05))
quantile(map_perc, seq(0, 1, by=0.05))

###
# SRP050272 - 83 samples, lncRNAs, some low map rates, high intergenic region counts
# ERP001908 - 62 miRNAs of tumors
# ERP001344 - 34 viral hepatitis agents
# ERP004592 - 23 miRNAs Hungtingtons
# SRP014675 - 43 miRNAs