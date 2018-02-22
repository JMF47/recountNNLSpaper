##### Contains the code to build all indices use in the analysis
cd /dcl01/leek/data/ta_poc/	
##########################################################################################
##### Full hg38 genome
##########################################################################################
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz -O /dcl01/leek/data/ta_poc/hisat2-indices/genome-index.tar.gz
cd /dcl01/leek/data/ta_poc/hisat2-indices/
tar -xvzf genome-index.tar.gz

##########################################################################################
##### Coding genes of GencodeV25 annotation
##########################################################################################

## HISAT2 ################################################################################
cd /dcl01/leek/data/ta_poc/GencodeV25
python ~/hisat2-2.0.5/hisat2_extract_exons.py gencodeV25.coding.gtf > gencodeV25.coding.exons
python ~/hisat2-2.0.5/hisat2_extract_splice_sites.py gencodeV25.coding.gtf > gencodeV25.coding.splice.sites
hisat2-build --exon gencodeV25.coding.exons --ss gencodeV25.coding.splice.sites \
	/dcl01/leek/data/ta_poc/hg38/hg38.fa \
	/dcl01/leek/data/ta_poc/hisat2-indices/hisat2-hg38-index

## Kallisto ##############################################################################
PATH=$PATH:/users/jmfu/kallisto_linux-v0.43.0/
index=kallisto/kallisto_index
kallisto index -i $index GencodeV25/gencodeV25.coding.fa

## Salmon ################################################################################
source /users/jmfu/tbb2018_20171205oss/bin/tbbvars.sh intel64 linux auto_tbbroot
cd ~/Salmon-0.8.2_linux_x86_64/bin/
out_dir=GencodeV25/salmon_index 
./salmon index -t GencodeV25/gencodeV25.coding.fa -i $out_dir --type quasi -k 31

## RSEM ##################################################################################
module load bowtie/1.1.1
PATH=$PATH:/users/jmfu/RSEM-1.3.0
cd GencodeV25
rsem-prepare-reference --gtf gencodeV25.coding.gtf --bowtie hg38.fa rsem-index/rsem

## sailfish ##############################################################################
PATH=$PATH:/users/jmfu/SailfishBeta-0.7.6_Linux-x86-64/bin
sailfish index -t GencodeV25/gencodeV25.coding.fa  -o  sailfish-index -k 31 -p 1

##########################################################################################
##### Coding genes [GencodeV25, chr1+chr14, simulations]
##########################################################################################

## Kallisto ##############################################################################
PATH=$PATH:/users/jmfu/kallisto_linux-v0.43.0/
kallisto index -i hg38_sim_kallisto_ind /dcl01/leek/data/ta_poc/geuvadis/GencodeV25/gencodeV25.sim.fa

## Salmon ################################################################################
source /users/jmfu/tbb2018_20171205oss/bin/tbbvars.sh intel64 linux auto_tbbroot
cd ~/Salmon-0.8.2_linux_x86_64/bin/
out_dir=/dcl01/leek/data/ta_poc/geuvadis/GencodeV25/salmon_index 
./salmon index -t /dcl01/leek/data/ta_poc/geuvadis/GencodeV25/sim.fa -i hg38_sim_salmon_ind --type quasi -k 31

## RSEM ##################################################################################
module load bowtie/1.1.1
PATH=$PATH:/users/jmfu/RSEM-1.3.0
rsem-prepare-reference --gtf gencodeV25.sim.gtf --bowtie hg38_sim.fasta rsem_sim_ind

