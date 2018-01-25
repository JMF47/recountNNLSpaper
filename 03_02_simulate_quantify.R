base_dir = "/dcl01/leek/data/ta_poc/geuvadis"
setwd(base_dir)
##### Quantify expression on simulated data in R
##### Currently set to quantify 75bp single end simulation
rl = 75
paired = 1
condition = paste0(rl, "_", paired)
samples = paste0("sample_0", 1)
outdir = paste0("simulation/", condition)
dir.create(outdir)
setwd(outdir)

##########################################################################################
# ### Rail-RNA
# if(paired==1){
# 	files = paste0("simulation/", condition, "/reads/", samples, ".fasta")
# 	man = cbind(files, 0, samples)
# 	rail_dir = paste0("simulation/", condition, "/rail")
# 	dir.create(rail_dir)
# 	write.table(man, sep="\t", quote=F, col.names=F, row.names=F, file=paste0(rail_dir, "/manifest.txt"))	
# }else{
# 	files = paste0("simulation/", condition, "/reads/", samples)
# 	man = cbind(paste0(files, "_1.fasta"), 0, paste0(files, "_2.fasta"), 0, samples)
# 	rail_dir = paste0("simulation/", condition, "/rail")
# 	dir.create(rail_dir)
# 	write.table(man, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, file=paste0(rail_dir, "/manifest.txt"))	
# }
# dir.create(paste0("simulation/", condition, "/recount/"))
# dir.create(paste0("simulation/", condition, "/recount/logs"))
# index="/dcl01/leek/data/HG38_rail/Homo_sapiens/UCSC/hg38/Sequence/"
# setwd(rail_dir)
# system2("rail-rna", paste0("go local -x ", index, "/BowtieIndex/genome ", index, "/Bowtie2Index/genome -m ", rail_dir, "/manifest.txt"))

setwd(base_dir)
library(recountNNLS)
bw = paste0('simulation/', condition, '/rail/rail-rna_out/coverage_bigwigs/sample_01.bw')
jx_file = paste0('simulation/', condition, '/rail/rail-rna_out/cross_sample_results/junctions.tsv.gz')
table = data.frame(project=rl, run=paired, bigwig_path=bw, rls=rl, paired_end=(paired==2))
pheno = processPheno(table)
rse_tx = recountNNLS(pheno, jx_file=jx_file, cores=20)
setwd(outdir)
save(rse_tx, file='recount/rse_tx.rda')

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

