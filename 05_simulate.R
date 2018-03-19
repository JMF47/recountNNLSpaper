##### Contains the code to simulate the full synthetic data (chr1 + chr 14)
##### Currently set to 75bp, single end simulation
##### Change values of rl and pe to simulate the other conditions in the study

##########################################################################################
### Calculate the simulation truth
##########################################################################################
tx_sim=import('gencodeV25.sim.gtf')
txs_sim = txs_sim[txs_sim$type=="transcript"]
sim_info = data.frame(transcript_id = txs_sim$transcript_id, gene_id = txs_sim$gene_id, reads = 0)

library(gtools)
### A wrapper function to do dirichlet draws on a gene by gene basis 
### to determine relative abundance of transcripts at a locus
simAbundance = function(gene, sim_info){
	sim_info_sub = sim_info[sim_info$gene_id == gene,]
	k = dim(sim_info_sub)[1]
	abundance = rdirichlet(1, rep(1/k, k))
	return(abundance)
}

### Calculate relative proportion of isoforms at each gene
set.seed(76) ## rl(75)+1
genes = unique(sim_info$gene_id)
abundances = sapply(genes, simAbundance, sim_info)
abundances = do.call(c, abundances)
sim_info$abundance = abundances

### Calculate total number of reads for a gene
tx_perGene = by(rep(1, dim(sim_info)[1]), droplevels(sim_info$gene_id), sum)
tx_perGene = data.frame(gene_id = names(tx_perGene), tx_count = as.numeric(tx_perGene))
set.seed(77) ## rl(75)+2
tx_perGene$geneReads = rnbinom(dim(tx_perGene)[1], 4, 0.01)
sim_info = merge(sim_info, tx_perGene, by="gene_id", sort=F)

### Calculate the number of reads for an isoform
### It is a product of the gene reads and distribution amongst isoforms
### This is passed on to polyester
sim_info$reads = round(sim_info$geneReads*sim_info$abundance)

##########################################################################################
### Simulate the reads using polyester
##########################################################################################
rl=75
pe=FALSE
outdir=paste0("simulation/reads/", rl, "_", pe*1+1)
	dir.create(outdir); setwd("outdir")

fasta = import("sim.fa", format="fasta")
	ind = which(sim_info$reads>0)
	fasta = fasta[ind]
sim_info = sim_info[ind,]
writeXStringSet(fasta, "tx.fasta")
set.seed(rl+1) ## rl(75)+1
simulate_experiment('tx.fasta', reads_per_transcript=sim_info$reads, 
	num_reps=c(2, 2), fold_changes=1, outdir="reads", paired=pe, readlen=rl)

base_dir = "/dcl01/leek/data/ta_poc/geuvadis"
setwd(base_dir)

##########################################################################################
### Quantify
##########################################################################################
rm(list=ls())
rl = 75
paired = 1
condition = paste0(rl, "_", paired)
samples = paste0("sample_0", 1)
outdir = paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/", condition)
dir.create(outdir)
setwd(outdir)

# ##########################################################################################
# ### Rail-RNA
# if(paired==1){
# 	files = paste0("simulation/", condition, "/reads/", samples, ".fasta")
# 	man = cbind(files, 0, samples)
# 	rail_dir = paste0("simulation/", condition, "/rail")
	#dir.create(rail_dir)
# 	write.table(man, sep="\t", quote=F, col.names=F, row.names=F, file=paste0(rail_dir, "/manifest.txt"))	
# }else{
# 	files = paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/", condition, "/reads/", samples)
# 	man = cbind(paste0(files, "_1.fasta"), 0, paste0(files, "_2.fasta"), 0, samples)
# 	rail_dir = paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/", condition, "/rail")
# 	dir.create(rail_dir)
# 	write.table(man, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, file=paste0(rail_dir, "/manifest.txt"))	
# }
# dir.create(paste0("simulation/", condition, "/recount/"))
# dir.create(paste0("simulation/", condition, "/recount/logs"))
# index="/dcl01/leek/data/HG38_rail/Homo_sapiens/UCSC/hg38/Sequence/"
# setwd(rail_dir)
# system2("rail-rna", paste0("go local -x ", index, "/BowtieIndex/genome ", index, "/Bowtie2Index/genome -m ", rail_dir, "/manifest.txt"))

### recountNNLS
library(recountNNLS); library(recountNNLSdata)
samps = "sample_01"
bw = paste0('/dcl01/leek/data/ta_poc/geuvadis/simulation/', condition, '/rail/rail-rna_out/coverage_bigwigs/', samps, '.bw')
jx_file = paste0('/dcl01/leek/data/ta_poc/geuvadis/simulation/', condition, '/rail/rail-rna_out/cross_sample_results/junctions.tsv.gz')
table = data.frame(project=rl, run=samps, bigwig_path=bw, rls=rl, paired_end=(paired==2))
pheno = processPheno(table)
rse_tx = recountNNLS(pheno, jx_file=jx_file, cores=20)
save(rse_tx, file=paste0("~/", rl, "_", paired, "_0.rda"))

##########################################################################################
### HISAT2-cufflinks (2.0.5-2.2.1)
setwd(outdir)
dir.create("hisat2")
dir.create("hisat2/aligned")
dir.create("hisat2/sorted")

for(s in samples){
	print(s); flush.console()
	outPath = paste0(outdir, "/hisat2/sorted/", s, ".cl.bam")
	outBam = paste0(outdir, "/hisat2/sorted/", s, ".cl")
	outSam = paste0(outdir, "/hisat2/aligned/", s, ".cl.sam")
	intermediateBam = paste0(outdir, "/hisat2/aligned/", s, "cl..bam")
	if(paired==1){
		readPath = paste0(outdir, "/reads/", s, ".fasta")
		system2("hisat2", paste0("--dta-cufflinks -f -x /dcl01/leek/data/ta_poc/hisat2-indices/hisat2-hg38-index -U ", readPath, " -S ", outSam))	
	}else{
		readPath = paste0(outdir, "/reads/", s)
		system2("hisat2", paste0("--dta-cufflinks -f -x /dcl01/leek/data/ta_poc/hisat2-indices/hisat2-hg38-index -1 ", readPath, "_1.fasta", " -2 ", readPath, "_2.fasta ", "-S ", outSam))
	}
	system2("samtools", paste0("view -b ", outSam, " > ", intermediateBam))
	system2("samtools", paste0("sort ", intermediateBam, " ", outBam))
}

cufflinks_dir = paste0(outdir, "/hisat2/cufflinks")
system2("mkdir", cufflinks_dir)
cufflinks_dir_default = paste0(cufflinks_dir, "/default")
system2("mkdir", cufflinks_dir_default)
annotation="/dcl01/leek/data/ta_poc/GencodeV25/gencodeV25.sim.gtf"

setwd(outdir)
for(s in samples){
	outBam = paste0(outdir, "/hisat2/sorted/", s, ".cl")
	workingDir = paste0(cufflinks_dir, "/", s)
	workingDir2 = paste0(cufflinks_dir_default, "/", s)
	system2("mkdir", workingDir)
	system2("mkdir", workingDir2)
	system2("cufflinks", paste0("-m 250 -s 25 --total-hits-norm --no-effective-length-correction --no-length-correction -o ",
		workingDir, " -G ", annotation, " ", outBam, ".bam"))
	system2("cufflinks", paste0("-m 250 -s 25 -o ",
		workingDir2, " -G ", annotation, " ", outBam, ".bam"))
}

##########################################################################################
### Kallisto (0.43.0)
kallisto_dir = paste0(outdir, "/kallisto")
system2("mkdir", kallisto_dir)
index="kallisto/kallisto_index"

for(s in samples){
	workingDir = paste0(kallisto_dir, "/", s)
	system2("mkdir", workingDir)
	if(paired==1){
		system2("kallisto", paste0("quant -i ", index, " -o ", workingDir, " --single -l 250 -s 25 ", outdir, "/reads/", s, ".fasta"))
	}else{
		system2("kallisto", paste0("quant -i ", index, " -o ", workingDir , " ", outdir, "/reads/", s, "_1.fasta ", outdir, "/reads/", s, "_2.fasta"))
	}
}

##########################################################################################
### Salmon (0.8.2)
index="GencodeV25/salmon_index "
salmon_dir = paste0(outdir, "/salmon")
dir.create(salmon_dir)
for(s in samples){
	if(paired==1){
		system2("/users/jmfu/Salmon-0.8.2_linux_x86_64/bin/salmon", 
			paste0("quant -i ", index, " -l A -r ", outdir, "/reads/", s, ".fasta -o ", salmon_dir, "/", s))	
	}
	else{
		system2("/users/jmfu/Salmon-0.8.2_linux_x86_64/bin/salmon", 
			paste0("quant -i ", index, " -l A -1 ", outdir, "/reads/", s, "_1.fasta -2 ", outdir, "/reads/", s, "_2.fasta -o ", salmon_dir, "/", s))	
	}
}

##########################################################################################
### RSEM (1.3.0)
rsem_dir = paste0(outdir, "/rsem")
dir.create(rsem_dir)
index = "rsem-index/rsem "

for(s in samples){
	if(paired==1){
		system2("rsem-calculate-expression", paste0("--no-qualities ", outdir, "/reads/", s, ".fasta ", index, rsem_dir, "/", s))
	}else{
		system2("rsem-calculate-expression", paste0("--no-qualities --paired-end ",  outdir, "/reads/", s, "_1.fasta ", outdir, "/reads/", s, "_2.fasta ", index, rsem_dir, "/", s))
	}
}

##########################################################################################
### Compile the information across methods
##########################################################################################
rm(list=ls())
library(Biostrings); library(GenomicFeatures); library(stringr); library(rtracklayer)
library(SummarizedExperiment)
rl = 37
paired = 1
condition = paste0(rl, "_", paired)
samples = paste0("sample_0", 1)
outdir = paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/", condition)
setwd(outdir)
load("reads/sim_counts_matrix.rda")
rownames(counts_matrix) = stringr::str_replace(stringr::str_extract(rownames(counts_matrix), "ENST.*gene"), " gene", "")
TxDb = makeTxDbFromGFF("/dcl01/leek/data/ta_poc/GencodeV25/gencodeV25.sim.gtf")
    truth = transcriptLengths(TxDb)
    truth$reads = 0
    mat = match(truth$tx_name, rownames(counts_matrix))
    truth$reads[!is.na(mat)] = counts_matrix[mat[!is.na(mat)],1]

### Load and parse the recountNNLS counts
# load("recount/rse_tx_simplese.rda")
load(paste0("~/", rl, "_", paired, "_0.rda"))
recountNNLS = assays(rse_tx)$counts[match(truth$tx_name, rownames(rse_tx)),]
	recountNNLS = recountNNLS
recountNNLSse = assays(rse_tx)$se[match(truth$tx_name, rownames(rse_tx)),]
recountNNLSscore = assays(rse_tx)$score[match(truth$tx_name, rownames(rse_tx)),]
recountNNLSdf = assays(rse_tx)$df[match(truth$tx_name, rownames(rse_tx)),]

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
save(info, recountNNLSse, recountNNLSscore, recountNNLSdf, truth, file=paste0("~/", rl, "_", paired, "_info.rda"))