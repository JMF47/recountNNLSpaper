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
.lnnls = function(data, matrix){
      mat = matrix(apply(matrix, 2, as.numeric), nrow=dim(matrix)[1])
            # penalty = 1+1*(data==0)
            # mat = mat*penalty
            # data = data*penalty
      mod=nnls::nnls(mat, data)
      return(mod)
}
.bootNNLS = function(nnls, matrix, alpha, boot=1000){
    res = nnls$residuals
    beta = nnls$x
    mat = matrix(apply(matrix, 2, as.numeric), nrow=dim(matrix)[1])
    bootBs = matrix(0, ncol=dim(mat)[2], nrow=boot)
    for(i in seq_len(boot)){
          set.seed(i)
          res_new = sample(res, length(res), replace=TRUE)
          data_new = mat %*% beta + res_new
          bootBs[i,] = nnls::nnls(mat, data_new)$x
    }
    bootl = apply(bootBs, 2, function(x) quantile(x, alpha/2))
    bootu = apply(bootBs, 2, function(x) quantile(x, 1-alpha/2))
    return(list(cbind(bootl, bootu)))
}
.bootNNLS2 = function(nnls, matrix, alpha, boot=1000){
	beta = nnls$x
	mat = matrix(apply(matrix, 2, as.numeric), nrow=dim(matrix)[1])
	bootBs = matrix(0, ncol=dim(mat)[2], nrow=boot)
    for(i in seq_len(boot)){
		set.seed(i)
		test = sapply(1:length(beta), function(i) rmultinom(1, beta[i], P[,i]))
		data_new = apply(test, 1, sum)
		bootBs[i,] = nnls::nnls(mat, data_new)$x
    }
    bootl = apply(bootBs, 2, function(x) quantile(x, alpha/2))
    bootu = apply(bootBs, 2, function(x) quantile(x, 1-alpha/2))
    return(list(cbind(bootl, bootu)))
}
.bootNNLS3 = function(nnls, matrix, alpha, res_samples, boot=1000){
	beta = nnls$x
	mat = matrix(apply(matrix, 2, as.numeric), nrow=dim(matrix)[1])
	bootBs = matrix(0, ncol=dim(mat)[2], nrow=boot)
    for(i in seq_len(boot)){
		set.seed(i)
		# res_new = res_samples[,sample(1:dim(res_samples)[2], 1)]
		res_new = res_samples[cbind(1:dim(res_samples)[1], sample(1:dim(res_samples)[2], dim(res_samples)[1], replace=T))]
		data_new = (nnls$fitted + res_new)
		bootBs[i,] = nnls::nnls(mat, data_new)$x
    }
    bootl = apply(bootBs, 2, function(x) quantile(x, alpha/2))
    bootu = apply(bootBs, 2, function(x) quantile(x, 1-alpha/2))
    return(list(cbind(bootl, bootu)))
}

data(matrix_37)
ems = matrix_37
locus = 6

boot_cov = NULL
anal_cov = NULL
tot_size = NULL
for(locus in 1:1000){
	message(locus)
	set.seed(locus)
	genes = as.character(g2l$gene[g2l$locus==locus])
	ems_sub = ems[match(genes, names(ems))]
	ems_sub = ems_sub[sapply(ems_sub, length)>0]

	if(length(ems_sub)>1)
	    P = .mergeP(ems_sub)
	if(length(ems_sub)==1)
	    P = ems_sub[[1]]	

	beta = round(gtools::rdirichlet(1, alpha=rep(1/dim(P)[2], dim(P)[2])) *rnbinom(1, 4, 0.01))
	counts_base = beta %*% t(P)
	counts = counts_base
		# counts = round(counts_base)
		# counts = rpois(length(counts), counts)
		# counts = counts + (rnorm(length(counts), 0, 0.01))
		counts = counts + abs(rnorm(length(counts), 0, 0.1))
		# counts[counts>0] = counts[counts>0] + abs(rnorm(sum(counts>0), 0, 0.1))
		# counts[counts>0] = counts[counts>0] + rpois(sum(counts>0), counts[counts>0])
		counts = as.numeric(counts)

	beta = as.numeric(beta)
	lm_info = .llm(counts, P)
    sebeta = lm_info[[2]]

    b = .lnnls(counts, P)
	# boot_mods = .bootNNLS(b, P)
	boot_mods = .bootNNLS2(b, P)
	bootl = sapply(boot_mods, function(x) x[,1])
	bootu = sapply(boot_mods, function(x) x[,2])
	boot_cov = c(boot_cov, sum(beta>=bootl & beta<=bootu))

	## analytic
	anall = b$x-1.69*sebeta
	analu = b$x+1.69*sebeta
	anal_cov = c(anal_cov, sum(beta>=anall & beta<=analu))

	tot_size = c(tot_size, dim(P)[2])
	message(sum(boot_cov)/sum(tot_size))
	message(sum(anal_cov)/sum(tot_size))
}

perf = NULL
for(i in 1:20){
	b = .lnnls(counts[,i], P)
	boot_mods = .bootNNLS(b, P, alpha=0.05)
	# boot_mods = .bootNNLS2(b, P, alpha=0.05)
	bootl = sapply(boot_mods, function(x) x[,1])
	bootu = sapply(boot_mods, function(x) x[,2])
	cbind(b$x, beta[,i], bootl, bootu)
	perf = c(perf, mean(beta[,i]>=bootl & beta[,i]<=bootu))
}

load("~/bootstrap_sample.rda")
mods = apply(counts, 2, .lnnls, P)
bs = sapply(mods, function(x) x$x)
fitted = apply(P, 2, as.numeric) %*% bs
res = counts - fitted

perf=NULL
for(i in 1:20){
	test = .bootNNLS3(mods[[i]], P, 0.05, res)	
	perf = c(perf, mean(beta[,i]>= test[[1]][,1] & beta[,i]<= test[[1]][,2]))
}





