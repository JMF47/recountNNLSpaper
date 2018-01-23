##### Contains the code to simulate the simulation data
##### Currently set to 75bp, single end simulation

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

