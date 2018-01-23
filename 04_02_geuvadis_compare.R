##### Quantify using recountNNLS
##### Then compile and compare quantification between methods

##########################################################################################
# recountNNLS
##########################################################################################
library(recountNNLS)
pheno = processPheno('ERP001942')
rse_tx = recountNNLS(pheno, cores=20)

##########################################################################################
# compile results
##########################################################################################
library(rtracklayer)
library(recountNNLS); library(SummarizedExperiment)
load(getRseTx('ERP001942'))

### Select a sample to compare
set.seed(463)
sample = sample(colnames(rse_tx), 1)
sample = sort(colnames(rse_tx))[1]
tx_info = rowData(rse_tx)
tx_list = rownames(rse_tx)

### compile information
rn = assays(rse_tx)$counts[,sample]
st = import(paste0('geuvadis/hisat2-stringtie/stringtie/', sample, '.st'), format="gtf")
cl = import(paste0('geuvadis/hisat2-cufflinks/cufflinks/', sample, '/transcripts.gtf'))
rsem = read.table(paste0('geuvadis/rsem/', sample, '.isoforms.results'), header=T)
kallisto = read.table(paste0('geuvadis/kallisto/', sample, '/abundance.tsv'), header=T)
salmon = read.table(paste0('geuvadis/salmon/', sample, '/quant.sf'), header=T)
out = rn
st = st[st$type=="transcript"]
	st_mat = match(tx_list, st$transcript_id)
	st=as.numeric(st$cov[st_mat])
	st = st*tx_info$tx_len/75
	out = cbind(out, st)
cl = cl[cl$type=="transcript"]
	cl_mat = match(tx_list, cl$transcript_id)
	cl_cov=as.numeric(cl$cov[cl_mat])
	cl_cts = cl_cov*tx_info$tx_len/75
	out = cbind(out, cl_cts)
rsem_mat = match(tx_list, rsem$transcript_id)
	out = cbind(out, rsem=rsem$expected_count[rsem_mat])
kallisto_mat = match(tx_list, kallisto$target_id)
	out = cbind(out, kl = kallisto$est_counts[kallisto_mat])
salmon_mat = match(tx_list, salmon$Name)
	out = cbind(out, sl = salmon$NumReads[salmon_mat])
	out = data.frame(out)

### Spearman's correlation
cors=cor(out, method="spearman", use="complete")

### Agreement between expressed transcripts
threshold=0
c21=sum((out$out>threshold)*(out$kl>threshold), na.rm=T)
c31=sum((out$out>threshold)*(out$st>threshold), na.rm=T)
c41=sum((out$out>threshold)*(out$cl>threshold), na.rm=T)
c51=sum((out$out>threshold)*(out$rsem>threshold), na.rm=T)
c61=sum((out$out>threshold)*(out$sl>threshold), na.rm=T)
	c32=sum((out$kl>threshold)*(out$st>threshold), na.rm=T)
	c42=sum((out$kl>threshold)*(out$cl>threshold), na.rm=T)
	c52=sum((out$kl>threshold)*(out$rsem>threshold), na.rm=T)
	c62=sum((out$kl>threshold)*(out$sl>threshold), na.rm=T)
		c43=sum((out$st>threshold)*(out$cl>threshold), na.rm=T)
		c53=sum((out$st>threshold)*(out$rsem>threshold), na.rm=T)
		c63=sum((out$st>threshold)*(out$sl>threshold), na.rm=T)
			c54=sum((out$cl>threshold)*(out$rsem>threshold), na.rm=T)
			c64=sum((out$cl>threshold)*(out$sl>threshold), na.rm=T)
				c65=sum((out$rsem>threshold)*(out$sl>threshold), na.rm=T)
