library(nnlasso)
library(selectiveInference)

.nnlasso=function(y, x){
	cv = cv.nnlasso(x, y, family="normal", plot=FALSE)
	las = nnlasso(x, y, family="normal")
	las = nnlasso(x, y, family="normal", path=FALSE, eps=1, lambda = cv$lambda, SE=TRUE, intercept=FALSE, normalize=FALSE)
	out = cbind(floor(las$coef[2,]), las$se)
	return(out)
	# ind = which.min(abs(las$lambdas - cv$lambda))
	# # coef = las$coef[ind,]
	# coef = las$coef[2,]
	# inference = fixedLassoInf(x, y, coef, cv$lambda, tol.beta=0.1, family="gaussian")
	# out = matrix(NA, ncol=3, nrow=dim(x)[2])
	# out[inference$vars,] = cbind(coef[inference$vars], inference$ci)
	# return(out)
}

P_binary = P>0
P_bin_sum = apply(P_binary, 1, sum)
P_weight = 1/P_bin_sum

out = list()
for(i in 1:100){
	out[[i]] = .nnlasso(Ys[,i]*P_weight, Pmat*P_weight)
}

coefs = sapply(out, function(x) x[5,1])
ses = sapply(out, function(x) x[5,2])
mean(abs(coefs-100)/ses>1)


coefs = sapply(out, function(x) x[5,1])
cil = sapply(out, function(x) x[5,2])
ciu = sapply(out, function(x) x[5,3])

#### Simulated Data
.mergeP = function(P_list){
      for(i in 1:length(P_list)){
            P_list[[i]]$rn = rownames(P_list[[i]])
      }
      P_out = P_list[[1]]
      for(i in 2:length(P_list)){
            P_out = merge(P_out, P_list[[i]], by="rn", all=TRUE)
      }
      P_out[is.na(P_out)] = 0
      rownames(P_out) = P_out$rn
      P_out$rn = NULL
      return(P_out)
}
.llm = function(data, matrix){
      # mat = matrix(apply(matrix, 2, as.numeric), nrow=dim(matrix)[1])
      mod=lm(data~.-1, data=matrix)
      mod_out = summary(mod)$coefficients
      hat = hatvalues(mod)
      return(list(mod_out[,1], mod_out[,2], hat))
}
.lnnls = function(data, matrix, boot=100){
      mat = matrix(apply(matrix, 2, as.numeric), nrow=dim(matrix)[1])
      mod=nnls::nnls(mat, data)
      return(mod)
}
.calculateReads = function(locus, g2l, ems, counts, power=1, verbose=FALSE){
	message(locus)
	b = NULL; Vb = NULL; colinear_info = NULL; boot0=NULL; bootsd=NULL
	boot5 = NULL; boot95=NULL
	genes = as.character(g2l$gene[g2l$locus==locus])
	ems_sub = ems[match(genes, names(ems))]
	ems_sub = ems_sub[sapply(ems_sub, length)>0]
	if(length(ems_sub)>1)
	    P = .mergeP(ems_sub)
	if(length(ems_sub)==1)
	    P = ems_sub[[1]]
	if(length(ems_sub)==0)
	    P = NULL
 	## Determine weighting of the uniqueness of features
 	mat = match(rownames(P), rownames(counts))
    P_size_orig = dim(P)[2]
	P = P[which(!is.na(mat)),,drop=FALSE]
	counts_sub = counts[mat[!is.na(mat)],, drop=FALSE]

	test = NULL
	### Multi transcript-location
	if(dim(P)[2]>1){
		lm_info = matrix(apply(counts_sub, 2, .llm, P))
		beta = sapply(lm_info, function(x) x[[1]])
		colinear = which(colnames(P) %in% rownames(beta)==FALSE)
		if(length(colinear)>0){ ## If there are colinear transcripts
			colinear_tx = colnames(P)[colinear]
			P = P[,-colinear]
			colinear_info = data.frame(locus_id=locus, transcript_id = colinear_tx,
			                         P_dim=P_size_orig, P_colin=length(colinear))
		}
		if(dim(P)[2]>1){ 
			test = .nnlasso(counts_sub[,1], apply(P, 2, as.numeric))
		}
	}
	return(test)
}

nnlass = mclapply(1:20, .calculateReads, g2l, matrix_list, counts, mc.cores=20)





