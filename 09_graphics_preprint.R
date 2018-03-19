#################################################################################
### Creating features
#################################################################################
par(xpd=NA)
layout(matrix(1:3, ncol=1), heights = c(2, 2, 3))
par(mar=c(2, 0.5, 2, 0.5))
plot(c(0, 900), c(-0.5,3), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="Annotated Features - Gene i", xaxs="i", yaxs="i")
rect(xleft=c(100, 500), xright=c(400, 800), ytop = 2.5, ybot=1.5)
segments(x0=400, x1=500, y0=2, y1=2)
text(x=250, y=2, labels = "Exon 1 - 300bp")
text(x=650, y=2, labels = "Exon 2 - 300bp")
rect(xleft=c(100), xright=300, ytop=1, ybot=0)
text(x=100, y=c(2, 0.5), labels = c("Tx 1", "Tx 2"), pos=2)
text(x=200, y=0.5, labels = "Exon 3 - 200bp")

polygon(x=c(425, 475, 450), y=c(-1, -1, -1.4), col=1)

plot(c(0, 900), c(-0.5,3), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="Disjoined Features - Gene i", xaxs="i", yaxs="i")
rect(xleft=c(100, 300, 500), xright=c(300, 400, 800), ytop = 2.5, ybot=1.5)
segments(x0=400, x1=500, y0=2, y1=2)
text(x=200, y=2, labels = "200bp")
text(x=350, y=2, labels = "100bp")
text(x=650, y=2, labels = "300bp")
rect(xleft=c(100), xright=300, ytop=1, ybot=0)
text(x=200, y=0.5, labels = "200bp")
text(x=100, y=c(2, 0.5), labels = c("Tx 1", "Tx 2"), pos=2)

polygon(x=c(425, 475, 450), y=c(-1, -1, -1.4), col=1)

plot(c(0, 900), c(-2.25,3), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="Final Features - Gene i", xaxs="i", yaxs="i")
rect(xleft=c(100, 200, 300, 500, 600, 700), xright=c(200, 300, 400, 600, 700, 800), ytop = 2.5, ybot=1.5)
segments(x0=400, x1=500, y0=2, y1=2)
text(c(150, 250, 350, 550, 650, 750), y=rep(2, 6), labels=rep("100bp", 6))
rect(xleft=c(100, 200), xright=c(200, 300), ytop=1, ybot=0)
text(c(150, 250), y=0.5, labels=rep("100bp", 2))
text(x=100, y=c(2, 0.5), labels = c("Tx 1", "Tx 2"), pos=2)

rect(xleft=c(100, 200, 300, 400, 500, 600, 700), xright=c(200, 300, 400, 500, 600, 700, 800), ytop = -.5, ybot=-1.5, lty=c(1, 1, 1, 2, 1, 1, 1), border=2)
text(x=100, y=c(-1), labels = c("Feature"), pos=2, col=2)
text(c(150, 250, 350, 450, 550, 650, 750), y=rep(-1, 7), labels=1:7, col=2)

#################################################################################
### Data generating process
#################################################################################
png("~/downloads/sim_robust_se/fig2.png", height=400, width=400)
layout(matrix(1:6, nrow=1), widths = c(1, 2, 2, 1, 0.5, 1))
par(mar=c(0, 0, 4, 0))
plot(c(0, 4), c(-0.5, 8), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
par(xpd=NA);text(x=2, y=8.5, labels="Gene i", cex=2, col=2); par(xpd=F)
text(x=2, y=seq(7, 1, by=-1), labels=paste0("Feature ", 1:7))
plot(c(0, 4), c(-0.5, 8), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="Feature \nProbabilities", xaxs="i", yaxs="i", bty="n")
text(x=1,y=seq(7,1,by=-1),labels=c(expression(X[1]^1), expression(X[2]^1),
                                   expression(X[3]^1), expression(X[4]^1), expression(X[5]^1), expression(X[6]^1), expression(X[7]^1)))
text(x=3,y=seq(7,1,by=-1),labels=c(expression(X[1]^2), expression(X[2]^2),
                                   expression(X[3]^2), expression(X[4]^2), expression(X[5]^2), expression(X[6]^2), expression(X[7]^2)))
par(xpd=NA)
text(x=1, y=8, labels="Tx 1")
text(x=3, y=8, labels="Tx 2")
text(x=1, y=7.6, labels=expression("(X"^1*")"))
text(x=3, y=7.6, labels=expression("(X"^2*")"))
lines(x=c(0,0), y=c(0.8, 7.2))
lines(x=c(0, 0.1), y=c(0.8, 0.8))
lines(x=c(0, 0.1), y=c(7.2, 7.2))
lines(x=c(4,4), y=c(0.8, 7.2))
lines(x=c(4, 3.9), y=c(0.8, 0.8))
lines(x=c(4, 3.9), y=c(7.2, 7.2))
par(xpd=F)
text(x=2, y=0, labels="X", font=2, cex=2)

plot(c(0, 4), c(-0.5, 8), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="True\nTranscript\nAbundances", xaxs="i", yaxs="i", bty="n")
text(x=0.5, y=c(0, 4), labels="X", cex=2, col=2)
text(x=1, y=4, labels="[", cex=2)
text(x=1.5, y=4, labels=expression(beta^1), cex=1.5)
text(x=2.5, y=4, labels=expression(beta^2), cex=1.5)
text(x=3, y=4, labels="]'", cex=2)
text(x=3.5, y=c(0, 4), labels="=", cex=2, col=2)
text(x=2, y=0, labels=expression(beta*"'"), cex=2)
text(x=1.5, y=4.5, labels="Tx 1")
text(x=2.5, y=4.5, labels="Tx 2")

plot(c(0, 4), c(-0.5, 8), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="Observed\nFeature\nCounts", xaxs="i", yaxs="i", bty="n")
text(x=2, y=seq(7, 1, by=-1), labels=c(expression(Y[1]), expression(Y[2]), expression(Y[3]), expression(Y[4]), expression(Y[5]), expression(Y[6]), expression(Y[7])))
lines(x=c(0.2,0.2), y=c(0.8, 7.2))
lines(x=c(0.2, 0.4), y=c(0.8, 0.8))
lines(x=c(0.2, 0.4), y=c(7.2, 7.2))
lines(x=c(3.8,3.8), y=c(0.8, 7.2))
lines(x=c(3.8, 3.6), y=c(0.8, 0.8))
lines(x=c(3.8, 3.6), y=c(7.2, 7.2))
text(x=2, y=0, labels=expression("Y"), cex=2)

plot(c(0, 1), c(-0.5, 8), xlab="", ylab="", xaxt="n", yaxt="n", type="n", xaxs="i", yaxs="i", bty="n")
text(x=0.5, y=c(0, 4), labels="+", cex=2, col=2)

plot(c(0, 4), c(-0.5, 8), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="Errors", xaxs="i", yaxs="i", bty="n")
text(x=2, y=seq(7, 1, by=-1), labels=c(expression(epsilon[1]), expression(epsilon[2]), expression(epsilon[3]), expression(epsilon[4]), expression(epsilon[5]), expression(epsilon[6]), expression(epsilon[7])))
lines(x=c(0.2,0.2), y=c(0.8, 7.2))
lines(x=c(0.2, 0.4), y=c(0.8, 0.8))
lines(x=c(0.2, 0.4), y=c(7.2, 7.2))
lines(x=c(3.8,3.8), y=c(0.8, 7.2))
lines(x=c(3.8, 3.6), y=c(0.8, 0.8))
lines(x=c(3.8, 3.6), y=c(7.2, 7.2))
text(x=2, y=0, labels=expression(epsilon), cex=2)
dev.off()

#################################################################################
### Performance in fully synthetic data
#################################################################################
rm(list=ls())
rls = c(37, 50, 75, 100, 150)
mae_complete = NULL
for(paired in 1:2){
	for(rl in rls){
		load(paste0("~/downloads/preprint_update/", rl, "_", paired, "_info.rda"))
		info = info
		truth = truth
		ses = recountNNLSse
		scores = recountNNLSscore
		df = recountNNLSdf
			df[is.na(df)] = 1
			df[df==0] = 1

		### Calculat error, as well as by score cat
		err = abs(info-truth$reads)
		score_cat = cut(scores, seq(0, 1, by=0.1))
		err_cats = apply(err, 2, function(x) by(abs(x), score_cat, mean, na.rm=T))
		mae_complete = rbind(mae_complete, apply(abs(err), 2, mean, na.rm=T))

		# png(paste0("~/downloads/preprint_update/graphics/", rl, "_", paired, ".png"), width=600, height=600)
		png(paste0("~/downloads/preprint_update/graphics/", rl, "_", paired, ".png"), width=900, height=300)
		cols = c("black", "red", "orange", "purple", "green")
		# par(mfrow=c(2,2))
		par(mfrow=c(1,3))
		par(ps = 12, cex = 1, cex.main = 1)
		par(mar=c(4.5, 4, 4, 2))
		hist(scores, breaks=seq(0, 1, by=0.1), xlab="Uniqueness", 
			main="Distribution of Uniqueness Scores", col=0)
		u = par()$usr; rect(u[1], u[3], u[2], u[4], col=rgb(0, 0, 0, 0.2), border=NA)
		text(x=u[2]-(u[2]-u[1])/20, y=u[4]-(u[4]-u[3])/20, label="(A)")
		hist(scores, breaks=seq(0, 1, by=0.1), add=T, col=0)

		### calculate CI coverage as a function of score cat
		stats = qt(0.975, df)
		cis = data.frame(cil=info[,1]-stats*ses, ciu=info[,1]+stats*ses)
		cis$truth = truth$reads; cis$bs = info[,1]; cis$scores = scores
		inds = which(cis$bs>0)
		cis = cis[inds,]
		score_cat_1 = cut(cis$scores, seq(0, 1, by=0.1))
		med_se = by(ses[inds], score_cat_1, median, na.rm=T)
		hit = (cis$cil<=cis$truth & cis$ciu>=cis$truth)
		means = by(hit, score_cat_1, mean)

		plot(err_cats[,1], ylab="Mean Absolute Error", xlab="", xaxt="n", ylim=c(0, max(err_cats)),
		main="Mean Absolute Error", pch=19, ty="b")
		u = par()$usr; rect(u[1], u[3], u[2], u[4], col=rgb(0, 0, 0, 0.2), border=NA)
		for(i in 2:5) points(err_cats[,i], col=cols[i], pch=19, ty="b")
		legend('topright', legend = c("rc", "kl", "cl", "rsem", "sl"), col=cols, fill=cols, ncol=2, cex=0.8)
		axis(side=1, at = 1:length(means), labels = names(means), las=2)
		text(x=mean(u[1:2]), y=u[3], labels="Uniqueness", pos=3)
		text(x=u[1]+(u[2]-u[1])/20, y=u[3]+(u[4]-u[3])/20, label="(B)")

		med_se = log(med_se+1, 10)
		plot(med_se, xaxt="n", ylab="log10 Median SE", xlab="", main="Median SE", pch=19, ylim=c(0, max(med_se)))
		u = par()$usr; rect(u[1], u[3], u[2], u[4], col=rgb(0, 0, 0, 0.2), border=NA)
		axis(side=1, at = 1:length(means), labels = names(means), las=2)
		text(x=mean(u[1:2]), y=u[3], labels="Uniqueness", pos=3)
		# text(x=u[2]-(u[2]-u[1])/20, y=u[4]-(u[4]-u[3])/20, label="(C)")
		text(x=u[1]+(u[2]-u[1])/20, y=u[3]+(u[4]-u[3])/20, label="(C)")
		dev.off()

		### 1-rep CIs
		# hit = (cis$cil<=cis$truth & cis$ciu>=cis$truth)
		# means1 = by(hit, score_cat_1, mean)
		# plot(means, xaxt="n", xlab="", ylab="Average 95% CI Coverage", pch=19,
		# ylim=c(0.8, 1), main="CI Coverage vs Uniqueness")
		# u = par()$usr; rect(u[1], u[3], u[2], u[4], col=rgb(0, 0, 0, 0.2), border=NA)
		# axis(side=1, at = 1:length(means), labels = names(means), las=2)
		# abline(h=0.95, col=2)
		# text(x=mean(u[1:2]), y=u[3], labels="Uniqueness", pos=3)
		# # text(x=u[1]+(u[2]-u[1])/20, y=u[3]+(u[4]-u[3])/20, label="(D)")
		# text(x=u[2]-(u[2]-u[1])/20, y=u[4]-(u[4]-u[3])/20, label="(D)")
		# dev.off()

  		### 100-rep CIs
  # 		if(paired==1){
  # 			load(paste0("/users/jackfu/downloads/CI_cal/CI_cal_countmat_", rl, ".rda"))
		# 	load(paste0("/users/jackfu/downloads/CI_cal/CI_", rl, ".rda"))	
		# }else{
		# 	load(paste0("/users/jackfu/downloads/CI_cal/CI_cal_countmat_", rl, "_2.rda"))
		# 	load(paste0("/users/jackfu/downloads/CI_cal/CI_", rl, "_2.rda"))	
		# }
		# rse_sub = rse_tx[match(rownames(count_mat), rownames(rse_tx)),]
		# bs = assays(rse_sub)$fragments
		# se = assays(rse_sub)$ses
		# score = assays(rse_sub)$scores
		# df = assays(rse_sub)$df
		# stats = qt(0.975, as.numeric(df[,1]))
		# cil = bs-stats*se
		# ciu = bs+stats*se
		# hit = (cil<=count_mat)*(ciu>=count_mat)

		# hits = apply(hit, 1, mean)
		# score_cat = cut(score[,1], seq(0, 1, by=0.1))
		# thresholds = seq(0.9, 1, by=0.01)
		# means = NULL
		# for(thresh in thresholds){
		# 	means = rbind(means,  by(hits, score_cat, function(x) mean(x<thresh, na.rm=T)))
		# }
		# plot(1-means[6,], xaxt="n", xlab="", ylab="Proportion of Transcripts", pch=19,
		# 	ylim=c(0, 1), main="Transcripts with Adequate 95% CIs")
		# u = par()$usr; rect(u[1], u[3], u[2], u[4], col=rgb(0, 0, 0, 0.2), border=NA)
		# axis(side=1, at = 1:length(levels(score_cat)), labels = levels(score_cat), las=2)
		# text(x=mean(u[1:2]), y=u[3], labels="Uniqueness", pos=3)
		# # text(x=u[2]-(u[2]-u[1])/20, y=u[4]-(u[4]-u[3])/20, label="(D)")
		# text(x=u[2]-(u[2]-u[1])/20, y=u[3]+(u[4]-u[3])/20, label="(D)")
		# dev.off()
	}
}

### 
png("~/downloads/preprint_update/graphics/performance_rl.png", height=400, width=800)
layout(matrix(c(1, 1, 2, 3, 4, 4), nrow=3, byrow=T), heights = c(0.5, 5, 0.5))
par(ps = 12, cex = 1.2, cex.main = 1.2)
par(mar=rep(0,4))
plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
	text(x=0.5, y=0, "Mean Relative Difference vs Read Length", cex=2, pos=3)
par(mar=c(2, 4, 2, 0))
plot(mae_complete[1:5,1]~rls, ylab="Mean Absolute Error", xlab="", xaxt="n",
	main="Single End", pch=19, ty="b", ylim=c(0, max(mae_complete[1:5,])))
	u = par()$usr; rect(u[1], u[3], u[2], u[4], col=rgb(0, 0, 0, 0.2), border=NA)
	for(i in 2:5) points(mae_complete[1:5,i]~rls, col=cols[i], pch=19, ty="b")
	text(x=(u[2] - (u[2]-u[1])/20), y=(u[4]-(u[4]-u[3])/20), pos=1, labels="(A)")
legend('bottomright', legend = c("rc", "kl", "cl"), col=cols, fill=cols, ncol=3, cex=0.8)
axis(side=1, at = rls, labels = rls)
par(mar=c(2, 0, 2, 4))
plot(mae_complete[6:10,1]~rls, ylab="", xlab="", xaxt="n", yaxt="n",
	main="Paired End", pch=19, ty="b", ylim=c(0, max(mae_complete[1:5,])))
	u = par()$usr; rect(u[1], u[3], u[2], u[4], col=rgb(0, 0, 0, 0.2), border=NA)
	for(i in 2:5) points(mae_complete[6:10,i]~rls, col=cols[i], pch=19, ty="b")
	text(x=(u[2] - (u[2]-u[1])/20), y=(u[4]-(u[4]-u[3])/20), pos=1, labels="(B)")
legend('bottomleft', legend = c("rsem", "sl"), col=cols, fill=cols[4:5], ncol=2, cex=0.8)
axis(side=1, at = rls, labels = rls)
axis(side=4)
par(mar=rep(0,4))
plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
	text(x=0.5, y=0.3, "Read Lengths", cex=1.4, pos=3)
dev.off()

#################################################################################
### Performance in RSEM
#################################################################################
load("~/downloads/preprint_update/rsem_based_0.rda")
err = abs(out-count_mat)
mae = apply(abs(err), 2, mean, na.rm=T)
	mae[5] = mae[4]; mae[4] = NA
mae_complete = rbind(mae_complete, mae)
err[is.na(err)] = 0
score_cat = cut(score, seq(0, 1, by=0.1))
err_cats = apply(err, 2, function(x) by(abs(x), score_cat, mean, na.rm=T))

# png("~/downloads/preprint_update/graphics/rsem_based.png", height=600, width=600)
png("~/downloads/preprint_update/graphics/rsem_based.png", height=300, width=900)
cols = c("black", "red", "orange", "green")
# par(mfrow=c(2,2))
par(mfrow=c(1,3))
par(ps = 12, cex = 1, cex.main = 1)
par(mar=c(4.5, 4, 4, 2))
hist(score, breaks=seq(0, 1, by=0.1), xlab="Uniqueness", 
	main="Distribution of Uniqueness Scores", col=0)
u = par()$usr; rect(u[1], u[3], u[2], u[4], col=rgb(0, 0, 0, 0.2), border=NA)
hist(score, breaks=seq(0, 1, by=0.1), add=T, col=0)
text(x=u[2]-(u[2]-u[1])/20, y=u[4]-(u[4]-u[3])/20, label="(A)")

### calculate errors and CI coverage as a function of uniquness
stats = qt(0.975, df)
cis = data.frame(cil=out[,1]-stats*se, ciu=out[,1]+stats*se)
cis$truth = count_mat; cis$bs = out[,1]; cis$scores = score
inds = which(cis$bs>0)
cis = cis[inds,]
score_cat_1 = cut(cis$scores, seq(0, 1, by=0.1))
med_se = by(se[inds], score_cat_1, median, na.rm=T)
hit = (cis[,1]<=cis[,3] & cis[,2]>=cis[,3])
means = by(as.numeric(hit), score_cat_1, mean, na.rm=T)

plot(err_cats[,1], ylab="Mean Absolute Error", xlab="", xaxt="n", ylim=c(0, max(err_cats, na.rm=T)),
main="Mean Absolute Error", pch=19, ty="b")
u = par()$usr; rect(u[1], u[3], u[2], u[4], col=rgb(0, 0, 0, 0.2), border=NA)
for(i in 2:4) points(err_cats[,i], col=cols[i], pch=19, ty="b")
legend('topright', legend = c("rc", "kl", "cl", "sl"), col=cols, fill=cols, ncol=2, cex=0.8)
axis(side=1, at = 1:length(means), labels = names(means), las=2)
text(x=mean(u[1:2]), y=u[3], labels="Uniqueness", pos=3)
# text(x=u[1]+(u[2]-u[1])/20, y=u[3]+(u[4]-u[3])/20, label="(B)")
text(x=u[2]-(u[2]-u[1])/20, y=u[4]-(u[4]-u[3])/20, label="(B)")

med_se = log(med_se+1, 10)
plot(med_se, xaxt="n", ylab="log10 Median SE", xlab="", main="Median SE", pch=19, ylim=c(0, max(med_se)))
u = par()$usr; rect(u[1], u[3], u[2], u[4], col=rgb(0, 0, 0, 0.2), border=NA)
axis(side=1, at = 1:length(means), labels = names(means), las=2)
text(x=mean(u[1:2]), y=u[3], labels="Uniqueness", pos=3)
# text(x=u[1]+(u[2]-u[1])/20, y=u[3]+(u[4]-u[3])/20, label="(C)")
text(x=u[2]-(u[2]-u[1])/20, y=u[4]-(u[4]-u[3])/20, label="(C)")

# plot(means, xaxt="n", xlab="", ylab="Average 95% CI coverage", pch=19,
# ylim=c(0.8, 1), main="CI Coverage vs Uniqueness")
# u = par()$usr; rect(u[1], u[3], u[2], u[4], col=rgb(0, 0, 0, 0.2), border=NA)
# axis(side=1, at = 1:length(means), labels = names(means), las=2)
# abline(h=0.95, col=2)
# text(x=mean(u[1:2]), y=u[3], labels="Uniqueness", pos=3)
# # text(x=u[1]+(u[2]-u[1])/20, y=u[3]+(u[4]-u[3])/20, label="(D)")
# text(x=u[2]-(u[2]-u[1])/20, y=u[4]-(u[4]-u[3])/20, label="(D)")
dev.off()

#################################################################################
### Performance table
#################################################################################
rownames(mae_complete) = c( paste0(c(37, 50, 75, 100, 150), " single"), 
	paste0(c(37, 50, 75, 100, 150), " paired"), "RSEM-based")
xtable::xtable(mae_complete)

#################################################################################
### ERR188410
#################################################################################
rm(list=ls())
load("~/downloads/preprint_update/ERR188410.rda")
tx_info = rowData(rse_tx)
gff = rtracklayer::import("~/downloads/gencodeV25.coding.dj.gff3")

plotTranscripts = function(gene, tx_list, gff, showUnique=T, fig_lab){
      anno = gff[gff$gene_id==gene]
      anno_exons = anno[anno$type=="exon"]
      col = countOverlaps(anno_exons, anno_exons)
      anno_exons$colz = (col==1)*2
      anno_exons$bordz = (col==1)*1+1
      txs = tx_list
      xleft = min(start(anno)); xright = max(end(anno))
      width = xright-xleft
      plot(c(xleft-width/10, xright+width/10), c(1, length(txs)+1), ty="n", ylab="", xlab="", yaxt="n", main="", xaxs="i", yaxs="i", xaxt="n", bty="n")
      for(i in 1:length(txs)){
            t = txs[i]
            sub = anno_exons[anno_exons$transcript_id==t]
            # lwds = ((end(sub)-start(sub))<20)*5+1
            lwds = 1
            rect(xleft=start(sub), xright=end(sub), ytop = i-0.2+0.5, ybot=i+0.2+0.5, col=sub$colz, border=sub$bordz, lwd=lwds)
      }
      junctions = GRangesList()
      for(i in 1:length(txs)){
            tmp = anno[anno$transcript_id==txs[i]]
            tmp = tmp[tmp$type=="exon"]
            tmp = sort(reduce(tmp))
            jx = GRanges(seqnames = seqnames(tmp)[1], IRanges(end(tmp)[-length(tmp)]+1, start(tmp)[-1]-1))
            junctions[[i]] = jx
      }

      base = unlist(junctions)
      base = unique(base)

      jx_counts = countOverlaps(base, unlist(junctions), type="equal")
      unique_jx = base[jx_counts==1]

      for(i in 1:length(txs)){
            ol = findOverlaps(unique_jx, junctions[[i]], type="equal")
            if(length(ol)>0){
                  tmp = unique_jx[queryHits(ol)]
                  for(j in 1:length(ol)){
                        lines(x=c(start(tmp[j]), end(tmp[j])), y=c(i+0.5, i+0.5), col=2, lty=2)
                  }
            }
      }
      coords = par()$usr
      par(xpd=NA)
      text(x = (coords[2]-coords[1])*0.2+coords[1], y=(coords[4]-coords[3])*-0.075+coords[3], labels=fig_lab, cex=2)
      par(xpd=F)
}

counts = assays(rse_tx)$counts[,2] 
se = assays(rse_tx)$se[,2]
score = assays(rse_tx)$score[,2]

png("~/downloads/preprint_update/graphics/se.png", width=1000, height=450, pointsize = 20)
layout(matrix(c(1, 1, 2, 3, 4, 4, 5, 6), ncol=2, byrow=T), heights = c(0.5, 2, 0.5, 2), widths=c(4, 1))
par(mar=rep(0, 4))
plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
text(x=0.5, y=0.8, pos=1, paste0("KLHL17"), cex=2)

gene = "ENSG00000187961.13"
inds = which(tx_info$gene_id == gene)
txs = tx_info$tx_name[inds]
ord = order(se[inds], decreasing=T)
ord2 = order(se[inds], decreasing=F)
inds = inds[ord]
plotTranscripts(gene, gff, tx_list=txs[ord2], fig_lab="")

plot(c(0,1), c(1, length(txs)+1), ty="n", ylab="", xlab="", yaxt="n", main="", xaxs="i", yaxs="i", xaxt="n")
text(signif(counts[inds], 3), x=0.2, y=seq(length(txs)+1, 2, by=-1)-0.5)
text(signif(se[inds], 3), x=0.5, y=seq(length(txs)+1, 2, by=-1)-0.5, col=2)
text(signif(score[inds], 3), x=0.8, y=seq(length(txs)+1, 2, by=-1)-0.5)
par(xpd=NA)
text(x=c(.2, .5, .8), y=7, pos=1, c(expression(beta), "SE", "Uniq"), cex=2)

plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
text(x=0.5, y=0.8, pos=1, paste0("G6PD"), cex=2)

gene = "ENSG00000160211.15"
inds = which(tx_info$gene_id == gene)
txs = tx_info$tx_name[inds]
ord = order(se[inds], decreasing=T)
ord2 = order(se[inds], decreasing=F)
inds = inds[ord]
plotTranscripts(gene, gff, tx_list=txs[ord2], fig_lab="")

plot(c(0,1), c(1, length(txs)+1), ty="n", ylab="", xlab="", yaxt="n", main="", xaxs="i", yaxs="i", xaxt="n")
text(signif(counts[inds], 3), x=0.2, y=seq(length(txs)+1, 2, by=-1)-0.5)
text(signif(se[inds], 3), x=0.5, y=seq(length(txs)+1, 2, by=-1)-0.5, col=2)
text(signif(score[inds], 3), x=0.8, y=seq(length(txs)+1, 2, by=-1)-0.5)
par(xpd=NA)
text(x=c(.2, .5, .8), y=14.2, pos=1, c(expression(beta), "SE", "Uniq"), cex=2)
dev.off()

