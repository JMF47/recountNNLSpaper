##### Quantify expression on simulated data in R and then collate then across methods
##### Currently set to quantif 75bp single end reads

##########################################################################################
### Run methods
##########################################################################################

rl = 75
paired = 1
condition = paste0(rl, "_", paired)
samples = paste0("sample_0", 1)
outdir = paste0("simulation/", condition)
dir.create(outdir)
setwd(outdir)

##########################################################################################
### Rail-RNA
if(paired==1){
	files = paste0("simulation/", condition, "/reads/", samples, ".fasta")
	man = cbind(files, 0, samples)
	rail_dir = paste0("simulation/", condition, "/rail")
	dir.create(rail_dir)
	write.table(man, sep="\t", quote=F, col.names=F, row.names=F, file=paste0(rail_dir, "/manifest.txt"))	
}else{
	files = paste0("simulation/", condition, "/reads/", samples)
	man = cbind(paste0(files, "_1.fasta"), 0, paste0(files, "_2.fasta"), 0, samples)
	rail_dir = paste0("simulation/", condition, "/rail")
	dir.create(rail_dir)
	write.table(man, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, file=paste0(rail_dir, "/manifest.txt"))	
}
dir.create(paste0("simulation/", condition, "/recount/"))
dir.create(paste0("simulation/", condition, "/recount/logs"))
index="/dcl01/leek/data/HG38_rail/Homo_sapiens/UCSC/hg38/Sequence/"
setwd(rail_dir)
system2("rail-rna", paste0("go local -x ", index, "/BowtieIndex/genome ", index, "/Bowtie2Index/genome -m ", rail_dir, "/manifest.txt"))

library(recountNNLS)
bw = paste0('simulation/', condition, '/rail/rail-rna_out/coverage_bigwigs/sample_01.bw')
jx_file = paste0('simulation/', condition, '/rail/rail-rna_out/cross_sample_results/junctions.tsv.gz')
table = data.frame(project=rl, run=paired, bigwig_path=bw, rls=rl, paired_end=(paired==2))
pheno = processPheno(table)
rse_tx = recountNNLS(pheno, jx_file=jx_file, cores=20)

save(rse_tx, file=paste0('simulation/', condition, '/recount/rse_tx.rda'))

##########################################################################################
### HISAT2-StringTie (2.05-1.3.3b)
setwd(outdir)
dir.create("hisat2")
dir.create("hisat2/aligned")
dir.create("hisat2/sorted")
for(s in samples){
	print(s); flush.console()
	outPath = paste0(outdir, "/hisat2/sorted/", s, ".ss.bam")
	outBam = paste0(outdir, "/hisat2/sorted/", s, ".ss")
	outSam = paste0(outdir, "/hisat2/aligned/", s, ".ss.sam")
	intermediateBam = paste0(outdir, "/hisat2/aligned/", s, "ss..bam")
	if(paired==1){
		readPath = paste0(outdir, "/reads/", s, ".fasta")
		system2("hisat2", paste0("--dta -f -x isat2-full-index/hisat2-ind-detailed -U ", readPath, " -S ", outSam))
	}else{
		readPath = paste0(outdir, "/reads/", s)
		system2("hisat2", paste0("--dta -f -x isat2-full-index/hisat2-ind-detailed -1 ", readPath, "_1.fasta", " -2 ", readPath, "_2.fasta ", "-S ", outSam))
	}
	system2("samtools", paste0("view -b ", outSam, " > ", intermediateBam))
	system2("samtools", paste0("sort ", intermediateBam, " ", outBam))
}

stringtie_dir = paste0(outdir, "/stringtie")
annotation="GencodeV25/gencodeV25.coding.gtf"
system2("mkdir", stringtie_dir)
system2("cd", stringtie_dir)

for(s in samples){
	workingDir = paste0(stringtie_dir, "/", s)
	system2("mkdir", workingDir)
	system2("stringtie", paste0("-o ", workingDir, "/cov-e.gtf -G ", annotation, " -e simulation/", condition, "/hisat2/sorted/", s, ".ss.bam") )
	system2("stringtie", paste0("-o ", workingDir, "/cov.gtf -G ", annotation, " simulation/", condition, "/hisat2/sorted/", s, ".ss.bam") )
}

##########################################################################################
### HISAT2-cufflinks (2.0.5)
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
		system2("hisat2", paste0("--dta-cufflinks -f -x hisat2-full-index/hisat2-ind-full -U ", readPath, " -S ", outSam))	
	}else{
		readPath = paste0(outdir, "/reads/", s)
		system2("hisat2", paste0("--dta-cufflinks -f -x hisat2-full-index/hisat2-ind-full -1 ", readPath, "_1.fasta", " -2 ", readPath, "_2.fasta ", "-S ", outSam))
	}
	system2("samtools", paste0("view -b ", outSam, " > ", intermediateBam))
	system2("samtools", paste0("sort ", intermediateBam, " ", outBam))
}

cufflinks_dir = paste0(outdir, "/hisat2/cufflinks")
system2("mkdir", cufflinks_dir)
cufflinks_dir_default = paste0(cufflinks_dir, "/default")
system2("mkdir", cufflinks_dir_default)
annotation="GencodeV25/gencodeV25.coding.gtf"

for(s in samples){
	workingDir = paste0(cufflinks_dir, "/", s)
	system2("mkdir", workingDir)
	system2("cufflinks", paste0("-m 250 -s 25 --total-hits-norm --no-effective-length-correction --no-length-correction -o ",
		workingDir, " -G ", annotation, " simulation/", condition, "/hisat2/sorted/", s, ".cl.bam"))
	system2("cufflinks", paste0("-m 250 -s 25 -o ",
		cufflinks_dir_default, " -G ", annotation, " simulation/", condition, "/hisat2/sorted/", s, ".cl.bam"))
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
### RSEM
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
### Compile quantification across results
##########################################################################################
setwd(outdir)
### Compile the true counts
load("reads/sim_counts_matrix.rda")
rownames(counts_matrix) = stringr::str_replace(stringr::str_extract(rownames(counts_matrix), "ENST.*gene"), " gene", "")
TxDb = makeTxDbFromGFF("~/gencodeV25.sim.gtf")
    truth = transcriptLengths(TxDb)
    truth$reads = 0
    mat = match(truth$tx_name, rownames(counts_matrix))
    truth$reads[!is.na(mat)] = counts_matrix[mat[!is.na(mat)],1]

recountNNLS = assays(rse_tx)$counts
	recountNNLS = recountNNLS[match(truth$transcript_id, rownames(recountNNLS)),]
	recountNNLS = recountNNLS/paired

kallisto = read.table(paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/", condition, "/kallisto/sample_01/abundance.tsv"), header=T)
kl_mat = match(truth$transcript_id, kallisto$target_id)
kl = kallisto$est_counts[kl_mat]

setwd(paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/", condition, "/stringtie"))
info = import("sample_01/cov-e.gtf")
info = info[info$type=="transcript"]
st_cov = data.frame(transcript_id = info$transcript_id, info$cov)
st_cov = st_cov[match(truth$transcript_id, st_cov$transcript_id),]
st = as.numeric(as.character(st_cov[,2]))*truth$tx_len/rl/paired

setwd(paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/", condition, "/salmon"))
tab = read.table("sample_01/quant.sf", header=T)
salmon = data.frame(transcript_id = tab[,1], reads = tab[,5])
salmon = salmon[match(truth$transcript_id, salmon$transcript_id),]
sl = salmon[,2]

setwd(paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/", condition, "/rsem"))
tab = read.table("sample_01.isoforms.results", header=T)
rsem = cbind(tab[,5])
rownames(rsem) = tab$transcript_id
rsem = rsem[match(truth$transcript_id, rownames(rsem)),]

tmp = import(paste0('/dcl01/leek/data/ta_poc/geuvadis/simulation/', condition, '/hisat2/cufflinks/sample_01/transcripts.gtf'))
cufflinks = tmp[tmp$type=="transcript"]
cufflinks_cov = data.frame(transcript_id = cufflinks$transcript_id, cufflinks$cov)
cufflinks_cov = cufflinks_cov[match(truth$transcript_id, cufflinks_cov$transcript_id),]
cl = as.numeric(as.character(cufflinks_cov[,2]))*truth$tx_len/rl/paired

info = data.frame(recountNNLS=recountNNLS, kl=kl, st=st, cl=cl, rsem=rsem, sl=sl)

