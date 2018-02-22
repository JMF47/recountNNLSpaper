rm(list=ls())
library(nnls)
library(recountNNLSdata)
data(matrix_75)
P = matrix_75[[1]]
Pmat = apply(P, 2, as.numeric)

### Normal distribution of errors on top of Y=Xbeta
sd=1
simulateY2 = function(beta, P){
      means = apply(t(P)*beta, 2, sum)
      Y = rnorm(length(means), means, sd=sd)
      # reads = sapply(1:length(beta), function(i) rnorm(dim(P)[1], P[,i]*beta[i], sd=0.1))
	# Y = apply(reads, 1, sum)
	return(Y)
}
sampleSD = function(values, boot=1000){
	set.seed(values[1])
	return(sd(sample(values, boot, replace=T)))
}
.lnnls = function(counts, P){
      # P_binary = P>0
      # P_bin_sum = apply(P_binary, 1, sum)
      # P_weight = 1/P_bin_sum
      # # P_weight = (P_bin_sum==1)*9+1
      # P = P*P_weight
      # counts = counts*P_weight
            # penalty = matrix(rep(apply(P, 2, mean)/40, dim(P)[1]), ncol = dim(P)[2], byrow=T)
            # P = P+penalty
      mat = matrix(apply(P, 2, as.numeric), nrow=dim(P)[1])
      mod=nnls::nnls(mat, counts)
      return(mod)
}
pw = function(P){
      P_binary = P>0
      P_bin_sum = apply(P_binary, 1, sum)
      P_weight = 1/P_bin_sum
      P = P*P_weight
      return(P)
}
r2 = function(P){
      r2 = sapply(1:dim(P)[2], function(i) summary(lm(Pmat[,i]~Pmat[,-i]-1))$'r.squared')
      return(r2)
}

set.seed(seed)
# beta = round(gtools::rdirichlet(1, alpha=rep(1/dim(P)[2], dim(P)[2])) *rnbinom(1, 4, 0.01))
# beta = rep(100, dim(P)[2])
# beta = c(100, 20, 50, 100, 150)
samp_size = 1000
beta = c(100, 0, 0, 0, 100)
results = NULL
results_cr = NULL
r2s = r2(Pmat)
n = dim(P)[1]
p = dim(P)[2]
fstat = qf(0.9, p, n-p)
threshold = p*fstat
Ys = matrix(ncol=samp_size, nrow=dim(P)[1])
for(j in 1:samp_size){
	set.seed(1+j)
	Ys[,j] = simulateY2(beta, P)
}
models = apply(Ys, 2, .lnnls, P)
bs = sapply(models, function(x) x$x)
      inc = apply(bs[2:4,], 2, sum)
      bs_correct = bs[,inc==0]
      sds_correct = sqrt(diag(cov(t(bs_correct)))[c(1, 5)])
      stat = (bs_correct[c(1, 5),] - 100)/sds_correct
      apply(stat, 1, function(x) mean(abs(x)<1.96))
ssRes = apply(sapply(models, function(x) x$residuals)^2, 2, sum)/(n-p)

covMat = cov(t(bs))
Pw = apply(pw(P), 2, as.numeric)
XTXinv = solve(t(Pw) %*% Pw)

## calculate betas, resample each row of beta to calculate sd, use sd to calculate CI
sds = apply(bs, 1, sampleSD)#*1/sqrt(1-r2s)
intl = bs - sds*1.96
intu = bs + sds*1.96
hitsl = matrix(rep(beta, samp_size), ncol=samp_size) >= intl
hitsu = matrix(rep(beta, samp_size), ncol=samp_size) <= intu

results = cbind(results, apply((hitsl*hitsu)[,1,drop=F], 1, mean))
results = apply((hitsl*hitsu), 1, mean)
stats = sapply(1:samp_size, function(i) rbind(bs[,i]-beta) %*% solve(covMat) %*% cbind(bs[,i]-beta))<threshold
stats = sapply(1:samp_size, function(i) rbind(bs[,i]-beta) %*% solve(XTXinv*ssRes[i]) %*% cbind(bs[,i]-beta))<threshold
k= which(beta==0)
stats_dropped = sapply(1:samp_size, function(i) 
      rbind(bs[-k,i]-beta[-k]) %*% solve(cov(t(bs[-k,]))) %*% cbind(bs[-k,i]-beta[-k]))<(qf(0.9, p-length(k), n-p)*(p-length(k)))


# results_cr = c(results_cr, rbind(bs[,1]-beta) %*% t(Pw) %*% Pw %*% cbind(bs[,1]-beta)<= threshold*sum(res[,1]^2)/(n-p))
test = rbind(apply(bs[c(1, 2, 4),], 2, sum), apply(bs[c(3, 5),], 2, sum))
      
      




