cd /dcl01/leek/data/ta_poc/
##### Script to quantify the geuvadis ERR188410 [excluding recountNNLS]
##### Paths to ERR188410

f=ERR188410
##########################################################################################
# HISAT2(2.0.5)-StringTie(1.3.3b)
##########################################################################################
## Index
cd GencodeV25
gffread -w hg38_coding.fa -g hg38.fa gencodeV25.coding.gtf
module load cufflinks/2.2.1 
python $hisat2_path/hisat2-2.0.5/hisat2_extract_exons.py gencodeV25.coding.gtf > gencodeV25.coding.exons
python $hisat2_path/hisat2-2.0.5/hisat2_extract_splice_sites.py gencodeV25.coding.gtf > gencodeV25.coding.splice.sites
hisat2-build --exon gencodeV25.coding.exons --ss gencodeV25.coding.splice.sites GencodeV25/hg38.fa hisat2-full-index/hisat2-ind-full

## Align + Quant
annotation=GencodeV25/gencodeV25.coding.gff3
cd geuvadis/
hisat2 --dta -p 20 -f -x hisat2-full-index/hisat2-ind-full -q -1 data/"$f"/"$f"_1.fastq.gz -2 data/"$f"/"$f"_2.fastq.gz -S hisat2-stringtie/"$f".sam
samtools view -b hisat2-stringtie/"$f".sam > hisat2-stringtie/"$f".bam
samtools sort hisat2-stringtie/"$f".bam -o hisat2-stringtie/sorted/"$f".bam -O 'bam' -T "stringtie_$f"
rm hisat2-stringtie/"$f".bam
rm hisat2-stringtie/"$f".sam
stringtie -o geuvadis/hisat2-stringtie/stringtie/"$f".st -e -G $annotation hisat2-stringtie/sorted/"$f".bam
rm hisat2-stringtie/sorted/"$f".bam

##########################################################################################
# HISAT2(2.0.5)-Cufflinks(2.2.1)
##########################################################################################
## Align + Quant
module load cufflinks/2.2.1
annotation=/dcl01/leek/data/ta_poc/GencodeV25/gencodeV25.coding.gff3
cd geuvadis
hisat2 --dta-cufflinks -p 20 -f -x /dcl01/leek/data/ta_poc/hisat2-indices/hisat2-hg38-index -q -1 data/"$f"/"$f"_1.fastq.gz -2 data/"$f"/"$f"_2.fastq.gz -S hisat2-cufflinks/"$f".sam
samtools view -b hisat2-cufflinks/"$f".sam > hisat2-cufflinks/"$f".bam
samtools sort hisat2-cufflinks/"$f".bam -o hisat2-cufflinks/sorted/"$f".bam -O 'bam' -T "cufflinks_$f"
rm hisat2-cufflinks/"$f".bam
rm hisat2-cufflinks/"$f".sam
cufflinks --total-hits-norm --no-effective-length-correction --no-length-correction \
-o /dcl01/leek/data/ta_poc/geuvadis/hisat2-cufflinks/cufflinks/"$f" \
-G /dcl01/leek/data/ta_poc/GencodeV25/gencodeV25.coding.gff3 \
/dcl01/leek/data/ta_poc/geuvadis/hisat2-cufflinks/sorted/"$f".bam
rm hisat2-cufflinks/sorted/"$f".bam

##########################################################################################
# RSEM(1.3.0 - bowtie 1.1.1)
##########################################################################################
module load bowtie/1.1.1
PATH=$PATH:/users/jmfu/RSEM-1.3.0
cd geuvadis
rsem-calculate-expression -p 10 --paired-end <(zcat data/"$f"/"$f"_1.fastq.gz) <(zcat data/"$f"/"$f"_2.fastq.gz) rsem-index/rsem rsem/"$f"

##########################################################################################
# Kallisto (0.43.0)
##########################################################################################
base=geuvadis/kallisto
fastq_base=geuvadis/data
index=geuvadis/kallisto/kallisto_index
cd geuvadis
kallisto quant -t 20 -i $index -o $base/$f $fastq_base/"$f"/"$f"_1.fastq.gz $fastq_base/"$f"/"$f"_2.fastq.gz

########################################################################
### Salmon (0.8.2)
########################################################################
source /users/jmfu/tbb2017_20161128oss/bin/tbbvars.sh intel64 linux auto_tbbroot
fastq_base=/dcl01/leek/data/ta_poc/geuvadis/data
index=/dcl01/leek/data/ta_poc/GencodeV25/hg38_coding_salmon_ind
cd ~/Salmon-0.8.2_linux_x86_64/bin/
./salmon quant -p 10 -i $index -l A -1 $fastq_base/$f/$f'_1.fastq.gz' -2 $fastq_base/$f/$f'_2.fastq.gz'  -o /dcl01/leek/data/ta_poc/geuvadis/salmon/$f

base_dir = "/dcl01/leek/data/ta_poc/"
setwd(base_dir)
##### Quantify using recountNNLS
##### Then compile and compare quantification between methods

##########################################################################################
# recountNNLS
##########################################################################################
library(recountNNLS)
pheno = processPheno('ERP001942')
	pheno = pheno[pheno$run %in% c("ERR188410", "ERR188412"),]
rse_tx = recountNNLS(pheno, cores=20)
save(rse_tx, file="~/ERR188410.rda")

##########################################################################################
# compile results
##########################################################################################
library(rtracklayer); library(recountNNLS); library(SummarizedExperiment)
load("~/ERR188410.rda")
sample = "ERR188410"
tx_info = rowData(rse_tx)
tx_list = rownames(rse_tx)

### compile information
rn = assays(rse_tx)$counts[,sample]
	scores = assays(rse_tx)$scores1[,sample]
cl = import(paste0('/dcl01/leek/data/ta_poc/geuvadis/hisat2-cufflinks/cufflinks/', sample, '/transcripts.gtf'))
rsem = read.table(paste0('/dcl01/leek/data/ta_poc/geuvadis/rsem/', sample, '.isoforms.results'), header=T)
kallisto = read.table(paste0('/dcl01/leek/data/ta_poc/geuvadis/kallisto/', sample, '/abundance.tsv'), header=T)
salmon = read.table(paste0('/dcl01/leek/data/ta_poc/geuvadis/salmon/', sample, '/quant.sf'), header=T)
out = rn
kallisto_mat = match(tx_list, kallisto$target_id)
	out = cbind(out, kl = kallisto$est_counts[kallisto_mat])
cl = cl[cl$type=="transcript"]
	cl_mat = match(tx_list, cl$transcript_id)
	cl_cov=as.numeric(cl$cov[cl_mat])
	cl_cts = cl_cov*tx_info$tx_len/75
	out = cbind(out, cl_cts)
rsem_mat = match(tx_list, rsem$transcript_id)
	out = cbind(out, rsem=rsem$expected_count[rsem_mat])
tx_list = stringr::str_replace(tx_list, "_PAR_Y", "")
salmon_mat = match(tx_list, salmon$Name)
	out = cbind(out, sl = salmon$NumReads[salmon_mat])
	out = data.frame(out)
out[is.na(out)] = 0

### Spearman's correlation
cors=cor(out, method="spearman")

### Agreement between expressed transcripts
out01 = out>0
t(out01) %*% out01

tx_info = rowData(rse_tx)
totals = apply(out, 2, sum)
fpkm = t(t(out/(tx_info$tx_len/1000))/(totals/1000000))

cor(out, method="spearman", use="complete")
cor(out[out[,1]>10,], method="spearman", use="complete")
out$test0 = out

cor(fpkm, method="spearman", use="complete")
cor(fpkm[fpkm[,1]>10,], method="spearman", use="complete")


