base_dir = "/dcl01/leek/data/ta_poc/geuvadis"
setwd(base_dir)
##### Collate quantified expression across methods
##### Currently set to quantif 75bp single end simulation
rl = 75
paired = 1
condition = paste0(rl, "_", paired)
samples = paste0("sample_0", 1)
outdir = paste0("simulation/", condition)
setwd(outdir)

##########################################################################################
### Compile quantification across results
##########################################################################################

### Compile the true counts
library(Biostrings); library(GenomicFeatures); library(stringr); library(rtracklayer)
load("reads/sim_counts_matrix.rda")
rownames(counts_matrix) = stringr::str_replace(stringr::str_extract(rownames(counts_matrix), "ENST.*gene"), " gene", "")
TxDb = makeTxDbFromGFF("/dcl01/leek/data/ta_poc/GencodeV25/gencodeV25.sim.gtf")
    truth = transcriptLengths(TxDb)
    truth$reads = 0
    mat = match(truth$tx_name, rownames(counts_matrix))
    truth$reads[!is.na(mat)] = counts_matrix[mat[!is.na(mat)],1]

### Load and parse the recountNNLS counts
# load("recount/rse_tx_simplese.rda")
load("recount/rse_tx_hybridse.rda")
recountNNLS = assays(rse_tx)$counts[match(truth$tx_name, rownames(rse_tx)),]
	recountNNLS = recountNNLS/paired
recountNNLSse = assays(rse_tx)$se[match(truth$tx_name, rownames(rse_tx)),]

### Load and parse the kallisto counts
kallisto = read.table("kallisto/sample_01/abundance.tsv", header=T)
kl_mat = match(truth$tx_name, kallisto$target_id)
kl = kallisto$est_counts[kl_mat]

### Load and parse the salmon counts
tab = read.table("salmon/sample_01/quant.sf", header=T)
salmon = data.frame(transcript_id = tab[,1], reads = tab[,5])
salmon = salmon[match(truth$tx_name, salmon$transcript_id),]
sl = salmon[,2]

### Load and parse the rsem counts
tab = read.table("rsem/sample_01.isoforms.results", header=T)
rsem = tab$expected_count
rsem = rsem[match(truth$tx_name, tab$transcript_id)]

### Load and parse the cufflinks counts
tmp = import('hisat2/cufflinks/sample_01/transcripts.gtf')
cufflinks = tmp[tmp$type=="transcript"]
cufflinks_cov = data.frame(transcript_id = cufflinks$transcript_id, cufflinks$cov)
cufflinks_cov = cufflinks_cov[match(truth$tx_name, cufflinks_cov$transcript_id),]
cl = as.numeric(as.character(cufflinks_cov[,2]))*truth$tx_len/rl/paired

### cbind all counts
info = data.frame(recountNNLS=recountNNLS, kl=kl, cl=cl, rsem=rsem, sl=sl)
# save(info, truth, recountNNLSse, file="results_simplese.rda")
save(info, truth, recountNNLSse, file="results_hybridse.rda")

##########################################################################################
### Compile metrics Across Scenarios
##########################################################################################
base_dir = "/dcl01/leek/data/ta_poc/geuvadis"
setwd(base_dir)
setwd("simulation")
table_rmse = NULL
table_mrd = NULL
ids = c("37_1", "50_1", "75_1", "100_1", "150_1", "37_2", "50_2", "75_2", "100_2", "150_2")
for(id in ids){
	load(paste0(id, "/results.rda"))
	info[is.na(info)] = 0
	err = info-truth$reads
		err[is.na(err)]=0
	rmse = sqrt(apply(err^2, 2, mean))
	table_rmse = rbind(table_rmse, rmse)
	rd = abs(info-truth$reads)/(info+truth$reads)*2
		rd[is.na(rd)] = 0
	# mrd = apply(rd*(truth$reads+1), 2, sum)/sum(truth$reads+1)
	mrd = apply(rd, 2, mean)
	table_mrd = rbind(table_mrd, mrd)
}
rownames(table_rmse)=ids
rownames(table_mrd)=ids
round(table_rmse, 1)
round(table_mrd, 3)
save(table_rmse, table_mrd, file="~/aggregate_results.rda")




