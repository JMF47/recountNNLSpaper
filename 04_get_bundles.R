##### Contains the code to determine bundles of genes to linear model simultaneously

rm(list=ls())
library(GenomicRanges); library(rtracklayer); library(stringr); library(Biostrings); library(Rsamtools)

### The exonic bins from our method for a read length
gff_ex = import('/dcl01/leek/data/ta_poc/GencodeV25/bins_150.bed')
### The probability matrices created previously for the same read length
load('/dcl01/leek/data/ta_poc/GencodeV25/matrix_150.rda')

### Calculation starts here
list_ex = rep( list(NULL), length(gff_ex))
list_jx = rep( list(NULL), length(gff_jx))

genes = names(mat)
for(i in 1:length(genes)){
	print(i); flush.console()
	gene = genes[i]
	tab = mat[[gene]]
	if(length(tab)>0){ # If the P is non-empty
		e_ind = str_replace(str_extract(rownames(tab), "e.*"), "e", "")
			e_ind = as.numeric(e_ind[!is.na(e_ind)])
		i_ind = str_replace(str_extract(rownames(tab), "i.*"), "i", "")
			i_ind = as.numeric(i_ind[!is.na(i_ind)])
		### Record the gene as having been observed in the corresponding exon and intron entries
		if(length(e_ind)>0){list_ex[e_ind] = lapply(list_ex[e_ind], function(x) c(x, gene))}
		if(length(i_ind)>0){list_jx[i_ind] = lapply(list_jx[i_ind], function(x) c(x, gene))}
	}
}

### Which genes need to be merged because of cross-mapping
num_genes_ex = sapply(list_ex, length)
num_genes_jx = sapply(list_jx, length)

list_ex1plus = list_ex[num_genes_ex>1]
	list_ex1p = unique(list_ex1plus)
list_jx1plus = list_jx[num_genes_jx>1]
	list_jx1p = unique(list_jx1plus)
list_1p = c(list_ex1p, list_jx1p)
	list_1p = unique(list_1p)

g2l = data.frame(locus = 1:length(genes), gene_id = genes)
counter = max(g2l$locus)+1
for(i in 1:length(list_1p)){
	message(i)
	list_sub = list_1p[[i]]
	ind = which(g2l$gene_id %in% list_sub)
	if(length(ind)>1){
		old_counter = g2l$locus[ind]
		if(length(unique(old_counter))>1){
			old_counter_ind = which(g2l$locus %in% old_counter)
			g2l$locus[old_counter_ind] = counter	
			counter = counter+1	
		}
	}
}
save(g2l_75, file="/dcl01/leek/data/ta_poc/GencodeV25/g2l_75.rda")
cts = by(rep(1, dim(g2l)[1]), g2l$locus, sum)
max(cts)
head(sort(cts, decreasing=T))

##########################################################################################
### 37 and 50bp where the networks create bundles of over 200 genes 
### 37bp is 4953
### 50bp is 1966 genes
##########################################################################################