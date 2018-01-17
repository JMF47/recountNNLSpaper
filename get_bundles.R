rm(list=ls())
library(GenomicRanges); library(rtracklayer); library(stringr); library(Biostrings); library(Rsamtools)
library(recountNNLSdata)
gff_ex = import('/dcl01/leek/data/ta_poc/GencodeV25/bins_50.bed')
load('/dcl01/leek/data/ta_poc/GencodeV25/matrix_50.rda')
mat = matrix_50

list_ex = rep( list(NULL), length(gff_ex))
list_jx = rep( list(NULL), length(gff_jx))

genes = names(mat)
for(i in 1:length(genes)){
	print(i); flush.console()
	gene = genes[i]
	tab = mat[[gene]]
	if(length(tab)>0){
		e_ind = str_replace(str_extract(rownames(tab), "e.*"), "e", "")
			e_ind = as.numeric(e_ind[!is.na(e_ind)])
		i_ind = str_replace(str_extract(rownames(tab), "i.*"), "i", "")
			i_ind = as.numeric(i_ind[!is.na(i_ind)])
		if(length(e_ind)>0){list_ex[e_ind] = lapply(list_ex[e_ind], function(x) c(x, gene))}
		if(length(i_ind)>0){list_jx[i_ind] = lapply(list_jx[i_ind], function(x) c(x, gene))}
	}
}

num_genes_ex = sapply(list_ex, length)
num_genes_jx = sapply(list_jx, length)

list_ex1 = list_ex[num_genes_ex==1]
list_ex1plus = list_ex[num_genes_ex>1]
list_jx1plus = list_jx[num_genes_jx>1]

gene_mat = data.frame(gene_id = genes, group=as.numeric(1:length(genes)))
for(i in 1:length(list_ex1plus)){
	message(i)
	working = list_ex1plus[[i]]
	groups = gene_mat$group[gene_mat$gene_id %in% working]
	gene_mat$group[gene_mat$gene_id %in% working] = rep(groups[1], length(groups))
}
for(i in 1:length(list_jx1plus)){
	message(i)
	working = list_jx1plus[[i]]
	groups = gene_mat$group[gene_mat$gene_id %in% working]
	gene_mat$group[gene_mat$gene_id %in% working] = rep(groups[1], length(groups))
}

gene_mat_50 = data.frame(locus=gene_mat[,2], gene = gene_mat[,1])
save(gene_mat_50, file="/dcl01/leek/data/ta_poc/GencodeV25/matrix_bundle_50.rda")

