##########################################################################################
### Preparing simulation truth (R)
##########################################################################################
rm(list=ls())
library(polyester); library(rtracklayer); library(Biostrings); library(stringr)
library(Biostrings); library(exomeCopy); library(GenomicFeatures); library(BSgenome.Hsapiens.UCSC.hg38)

### Obtain transcript sequences
hg38 = BSgenome.Hsapiens.UCSC.hg38
TxDb = makeTxDbFromGFF("/dcl01/leek/data/ta_poc/GencodeV25/gencodeV25.coding.gff3")
	TxL = transcriptLengths(TxDb)
	TxSeq = extractTranscriptSeqs(hg38, TxDb)

sim_info = data.frame(transcript_id = TxL$tx_name, gene_id = TxL$gene_id)
library(gtools)
simp = function(gene, sim_info){
	message(which(genes==gene))
	sim_info_sub = sim_info[sim_info$gene_id == gene,,drop=FALSE]
	k = dim(sim_info_sub)[1]
	if(k>1){
		abund = rdirichlet(1, rep(1/k, k))
	}else{
		abund=1
	}	
	return(abund)
}
seed=76
set.seed(seed)
	genes = unique(sim_info$gene_id)
	abundances = mclapply(genes, simp, sim_info, mc.cores=20)
	ps = do.call(c, abundances)
set.seed(seed+1)
	geneReads = rnbinom(length(genes), 4, 0.01)
	geneExpressed = rbinom(length(genes), 1, 0.82)
	gene_info = data.frame(gene_id = genes, geneReads, geneExpressed)
	sim_info = merge(sim_info, gene_info, by="gene_id", sort=F)
	sim_info$tx_p = ps
		test = by(sim_info$tx_p, sim_info$gene_id, sum)

gene_abundances = nbinom

sim_info$abundance = abundances

tx_perGene = by(rep(1, dim(sim_info)[1]), droplevels(sim_info$gene_id), sum)
tx_perGene = data.frame(gene_id = names(tx_perGene), tx_count = as.numeric(tx_perGene))
set.seed(seed+1)
tx_perGene$geneReads = rnbinom(dim(tx_perGene)[1], 4, 0.01)
sim_info = merge(sim_info, tx_perGene, by="gene_id", sort=F)
sim_info$reads = round(sim_info$geneReads*sim_info$abundance)
write.table(sim_info, row.names=F, file="sim_info.txt")


