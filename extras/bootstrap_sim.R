rm(list=ls())
library(nnls)
library(recountNNLSdata)
data(matrix_75)
P = matrix_75[[1]]
Pmat = apply(P, 2, as.numeric)

### Multinomial basis for simulation
simulateY2 = function(beta, P){
	means = apply(t(P)*beta, 2, sum)
	Y = rnorm(length(means), means, sd=1)
	# reads = sapply(1:length(beta), function(i) rnorm(dim(P)[1], P[,i]*beta[i], sd=0.1))
	# Y = apply(reads, 1, sum)
	return(Y)
}
sampleSD = function(values, boot=1000){
	set.seed(values[1])
	return(sd(sample(values[values>0], boot, replace=T)))
}
.lnnls = function(data, matrix){
      mat = matrix(apply(matrix, 2, as.numeric), nrow=dim(matrix)[1])
      mod=nnls::nnls(mat, data)
      return(mod)
}

### Simulate the betas 
beta = c(100, 100, 0, 100, 100)

### Simulate a dataset of size samp_size, with eacn sample sharing the common beta
set.seed(1337)
iter = 100
results = NULL
results2 = NULL
called = NULL
for(i in 1:iter){
	samp_size = 20
	Ys = matrix(ncol=samp_size, nrow=dim(P)[1])
	for(j in 1:samp_size){
		set.seed(i+j)
		Ys[,j] = simulateY2(beta, P)
	}
	models = apply(Ys, 2, .lnnls, P)
	bs = sapply(models, function(x) x$x)

	### Bootstrap #1
	## calculate betas
	## resample each row of beta to calculate sd
	## use sd to calculate CI
	sds = apply(bs, 1, sampleSD)
	intl = bs - sds*1.69
	intu = bs + sds*1.69
	hitsl = (matrix(rep(beta, samp_size), ncol=samp_size) >= intl)[,1]
	hitsu = (matrix(rep(beta, samp_size), ncol=samp_size) <= intu)[,1]
	called = cbind(called, bs[,1]>0)
	results = cbind(results, hitsl*hitsu)
}
apply(results*called, 1, sum)/apply(called, 1, sum)

boxplot(t(results), xlab="Transcript", ylab="Coverage")
abline(h=0.9, col=2)


	# ### Bootstrap #2
	# ## residual sampling
	# ## residuals are sample across features within a sample
	# .bootNNLS = function(nnls, matrix, alpha, boot=1000){
	#     res = nnls$residuals
	#     beta = nnls$x
	#     mat = matrix(apply(matrix, 2, as.numeric), nrow=dim(matrix)[1])
	#     bootBs = matrix(0, ncol=dim(mat)[2], nrow=boot)
	#     for(i in seq_len(boot)){
	#           set.seed(i)
	#           res_new = sample(res, length(res), replace=TRUE)
	#           data_new = mat %*% beta + res_new
	#           bootBs[i,] = nnls::nnls(mat, data_new)$x
	#     }
	#     bootl = apply(bootBs, 2, function(x) quantile(x, alpha/2))
	#     bootu = apply(bootBs, 2, function(x) quantile(x, 1-alpha/2))
	#     return(list(cbind(bootl, bootu)))
	# }
	# boots = sapply(models, .bootNNLS, P, alpha = 0.1)
	# intl = sapply(boots, function(x) x[,1])
	# intu = sapply(boots, function(x) x[,2])
	# hitsl = matrix(rep(beta, samp_size), ncol=samp_size) >= intl
	# hitsu = matrix(rep(beta, samp_size), ncol=samp_size) <= intu
	# results2 = cbind(results2, apply(hitsl*hitsu, 1, mean))
}
