rm(list=ls())
library(nnls)
library(recountNNLSdata)
data(matrix_75)
P = matrix_75[[1]]
Pmat = apply(P, 2, as.numeric)

sd=2
simulateY2 = function(beta, P){
    # reads = sapply(1:length(beta), function(i) rnorm(dim(P)[1], P[,i]*beta[i], sd=1))
    # Y = apply(reads, 1, sum)
    means = apply(t(P)*beta, 2, sum)
	Y = rnorm(length(means), means, sd=sd)
	return(Y)
}
sampleSD = function(values, boot=1000){
	set.seed(values[1])
	return(sd(sample(values, boot, replace=T)))
}
.lnnls = function(counts, P){
      mat = matrix(apply(P, 2, as.numeric), nrow=dim(P)[1])
      mod=nnls::nnls(mat, counts)
      return(mod)
}

set.seed(1337)
beta = c(0, 70, 0, 0, 70)
samp_size = 500
Ys = matrix(ncol=samp_size, nrow=dim(P)[1])
for(j in 1:samp_size){
	set.seed(1+j)
	Ys[,j] = simulateY2(beta, P)
}

# models = apply(Ys, 2, .lnnls, P)
# bs = sapply(models, function(x) x$x)
# res = sapply(models, function(x) x$residuals)
# covMat = cov(t(bs[bhat>0,]))
# Sigma = cov(t(res))
Sigma = diag(dim(Ys)[1])

alpha = 0.05
funcL = function(mu){
	(pnorm((x-mu)/s) - pnorm((a-mu)/s))/(pnorm((b-mu)/s) - pnorm((a-mu)/s)) - (1-alpha/2)
} 
funcU = function(mu){
	(pnorm((x-mu)/s) - pnorm((a-mu)/s))/(pnorm((b-mu)/s) - pnorm((a-mu)/s)) - (alpha/2)
} 

cis = NULL
pis = NULL
results = NULL
for(i in 1:samp_size){
	y = Ys[,i]
	mod = .lnnls(Ys[,i], P)
	bhat = mod$x
	if(sum(bhat>0)>1){
		Xs = Pmat[,bhat>0]
		Xcs = Pmat[,bhat==0]
		Xss = solve(t(Xs) %*% Xs) %*% t(Xs)

		### eta
		etas = t(Xss)
		for(j in seq_len(sum(bhat>0))){
			eta = etas[,j,drop=FALSE]
			x = (t(etas)%*%y)[j]
			a = 0
			b = sum(bhat)
			mu = beta[which(bhat>0)[j]]
			s = sqrt(sd^2*sum(eta^2))
			results = c(results, (pnorm((x-mu)/s) - pnorm((a-mu)/s))/(pnorm((b-mu)/s) - pnorm((a-mu)/s)))

			cil=uniroot(funcL, lower=0, upper=1, extendInt="yes")
			ciu=uniroot(funcU, lower=0, upper=1, extendInt="yes")
			cis = rbind(cis, c(x, mu, cil$root, ciu$root))
			pis = rbind(pis, c(x/b, mu/sum(beta), cil$root/b, ciu$root/b))
		}
	}
}
cis[cis[,3]<0,3]=0
cis_ind = (cis[,2]>=cis[,3] & cis[,2]<=cis[,4])

pis[pis[,3]<0,3]=0
pis_ind = (pis[,2]>=pis[,3] & pis[,2]<=pis[,4])




### Constraint Ay<=0
# A = rbind(		t(Xs) %*% (diag(dim(P)[1])- Xs %*% Xss),
# 				-t(Xs) %*% (diag(dim(P)[1])- Xs %*% Xss),
# 				-t(Xcs) %*% (diag(dim(P)[1])- Xs %*% Xss)
# 		)
# A = rbind(		t(Xs) %*% (diag(dim(P)[1])- Xs %*% Xss),
# 				-t(Xs) %*% (diag(dim(P)[1])- Xs %*% Xss)
# 		)
# A = rbind(-t(Xs) %*% (diag(dim(P)[1])- Xs %*% Xss))

# Adenom = as.numeric(t(eta) %*% Sigma %*% eta)
# alpha = (A %*% Sigma %*% eta)/Adenom

# Ay = A %*% y
# 	# Ay[1:4] = 0
# V = (-Ay + alpha * as.numeric(t(eta) %*% y))/alpha
# Vneg = max(V[alpha<0])
# Vpos = min(V[alpha>0])
# V0 = -Ay[alpha==0]
# a = max(0, Vpos)
# # a = Vpos
# b = min(Vneg, sum(bhat))