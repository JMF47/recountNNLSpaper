rm(list=ls())
load("~/downloads/sim_robust_se/truth.rda")
# load("~/downloads/150_1_rsem_HC3.rda")
# load("~/downloads/sim_robust_se/rsem_sim_75_2_1.rda")
load("~/downloads/sim_robust_se/rsem_sim_75_2_1_HC4.rda")
# load("~/downloads/sim_robust_se/test.rda")
# load("~/downloads/75_1_0_rsem.rda")
# 	truth = TxL
# 	ind = match(truth$tx_name, rownames(rse_tx))
# 	rse_sub = rse_tx[ind,]
# 	output1 = cbind(assays(rse_sub)$counts, assays(rse_sub)$se, assays(rse_sub)$scores)
# load("~/downloads/150_0_rsem_HC4m.rda")
# sl = read.table("~/downloads/quant.sf", header=TRUE)

ind = match(rownames(count_mat), rownames(rse_tx))
rse_sub = rse_tx[ind,]
output = cbind(assays(rse_sub)$counts, assays(rse_sub)$se, assays(rse_sub)$scores)
bs = output[,1]
ses = output[,2]
scores = output[,3]
count_mat = count_mat*2
sl = output

# cor(output[,1], output1[,1], method="spearman", use="complete")
# threshold = 100
# select = which(output[,1]>threshold & output1[,1]>threshold)
# cor(output[select,1], output1[select,1], method="spearman", use="complete")
# plot(log(output[select,1]+1, 2)~log(output1[select, 1]+1, 2), pch=19, col=rgb(0,0,0,0.1), cex=0.2)
# vw = cbind(output[,1], output1[,1], truth$counts)

cor(output[,1], count_mat, method="spearman", use="complete")
select = which(output[,1]>1)
cor(output[select,1], count_mat[select], method="spearman", use="complete")

### Checking calibration of confidence intervals
covs = seq(0.5, 0.99, by=0.01)
stats = qnorm(covs)
obs_covs = NULL
expressed = which(output[,1]>0)
for(i in 1:length(stats)){
cis = cbind(output[,1]-stats[i]*output[,2], output[,1]+stats[i]*output[,2])
obs_covs = c(obs_covs, sum(count_mat>=cis[,1] & count_mat<=cis[,2] & output[,1]>0, na.rm=T)/length(expressed))
}
covs = 1-(1-covs)*2

info = cbind(output[,1], output[,1])

err = info-as.numeric(count_mat)
apply(abs(err), 2, mean, na.rm=T)
sqrt(apply(err^2, 2, mean, na.rm=T))
# err = output[,1]-truth$counts
# mean(abs(err), na.rm=T)
# sqrt(mean(err^2, na.rm=T))

threshold = 50
cols = c("black", "red", "orange", "purple", "green")
par(mfrow=c(2,2))
par(mar=c(4.5, 4, 4, 4))
hist(scores[output[,1]>0], breaks=seq(0, 1, by=0.1), xlab="Uniqueness", main="Distribution of Uniqueness Scores",
col=0)
u = par()$usr; rect(u[1], u[3], u[2], u[4], col=rgb(0, 0, 0, 0.2), border=NA)
hist(scores[output[,1]>0], breaks=seq(0, 1, by=0.1), add=T, col=0)
# hist(scores[output[,1]>threshold], breaks=seq(0, 1, by=0.1), col=1, add=T)

# plot(lkps_err[,1]~lkps_med, ylab="Mean Absolute Error", xlab="Uniqueness", 
#       main="Mean absolute error vs Uniqueness", ylim=c(0, max(lkps_err)), pch=19, ty="b")
#       u = par()$usr; rect(u[1], u[3], u[2], u[4], col=rgb(0, 0, 0, 0.2), border=NA)
#       for(i in 2:5) points(lkps_err[,i]~lkps_med, col=cols[i], pch=19, ty="b")
#       legend('topright', legend = c("rc", "kl", "cl", "rsem", "sl"), col=cols, fill=cols, ncol=2, cex=0.8)

cis = data.frame(cil=output[,1]-1.96*output[,2], ciu=output[,1]+1.96*output[,2])
cis$truth = count_mat
cis$bs = output[,1]
cis$scores = scores
cis$fpkm = fpkm
cis = cis[complete.cases(cis),]
score_cat = cut(cis$scores, seq(0, 1, by=0.1))
med_se = by((cis$ciu-cis$cil)/1.96, score_cat, median)

	mean(cis$truth>0 & cis$cil>0, na.rm=T)
	mean(cis$truth>0)

err_sub = err[match(rownames(cis), rownames(err)),]
abs_errs = apply(err_sub, 2, function(x) by(abs(x), score_cat, mean, na.rm=T))

inds = which(cis$bs>0)
cis = cis[inds,]
hit = as.numeric(cis$cil<=cis$truth & cis$ciu>=cis$truth)
score_cat = score_cat[inds]
means = by(hit, score_cat, mean, na.rm=T)

	# power = (cis$truth>0 & cis$cil>0)
	# powers = by(power, score_cat, mean, na.rm=T)

inds = which(cis$bs>threshold)
cis = cis[inds,]
score_cat = cut(cis$scores, seq(0, 1, by=0.1))
med_se1 = by((cis$ciu-cis$cil)/1.96, score_cat, median)
# err_sub1 = err[match(rownames(cis), rownames(err)),]
# abs_errs1 = apply(err_sub1, 2, function(x) by(abs(x), score_cat, mean))

plot(abs_errs[,1], ylab="Mean Absolute Error", xlab="Uniqueness", xaxt="n",
main="Mean absolute error vs Uniqueness", ylim=c(0, max(abs_errs)), pch=19, ty="b")
u = par()$usr; rect(u[1], u[3], u[2], u[4], col=rgb(0, 0, 0, 0.2), border=NA)
# for(i in 2:2) points(abs_errs[,i], col=cols[i], pch=19, ty="b")
# legend('topright', legend = c("rc", "kl", "cl", "rsem", "sl"), col=cols, fill=cols, ncol=2, cex=0.8)
# axis(side=1, at = 1:length(means), labels = names(means), las=2)


med_se = log(med_se+1, 10)
med_se1 = log(med_se1+1, 10)
plot(med_se, xaxt="n", ylab="log10 Median SE", xlab="", main="Median SE vs Uniqueness", ylim=c(0, max(c(med_se, med_se1))))
u = par()$usr; rect(u[1], u[3], u[2], u[4], col=rgb(0, 0, 0, 0.2), border=NA)
points(med_se, col=0, pch=19)
points(med_se)
# points(med_se1, pch=19)
# legend('bottomleft', legend = c(">0 reads", "> reads"), fill = 0:1)
axis(side=1, at = 1:length(means), labels = names(means), las=2)

hit = as.numeric((cis$cil<=cis$truth & cis$ciu>=cis$truth))
means1 = by(hit, score_cat, mean)
plot(means, xaxt="n", xlab="", ylab="Average 95% CI coverage",
ylim=c(0.8, 1), main="CI coverage vs Uniqueness")
u = par()$usr; rect(u[1], u[3], u[2], u[4], col=rgb(0, 0, 0, 0.2), border=NA)
points(means, col=0, pch=19)
points(means)
axis(side=1, at = 1:length(means), labels = names(means), las=2)
# points(means1, pch=19)
abline(h=0.95, col=2)
# legend('bottomleft', legend = c(">0 reads", "> reads"), fill = 0:1)




