################################################################################################################################################
### fig1.png ###################################################################################################################################
################################################################################################################################################
png("~/downloads/recount_graphics/fig1.png", width=700, height=700, units="px", pointsize=24)
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
dev.off()

################################################################################################################################################
### fig2.png ###################################################################################################################################
################################################################################################################################################
png("~/downloads/recount_graphics/fig2.png", width=700, height=500, units="px", pointsize=24)
layout(matrix(1:4, nrow=1), widths = c(1, 2, 2, 1))
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
dev.off()

################################################################################################################################################
### Sim results
################################################################################################################################################
cs = 0.1
ids = c("37_1", "37_2", "50_1", "50_2", "75_1", "75_2", "100_1", "100_2", "150_1", "150_2")
subset = FALSE
for(id in ids){
      setwd(paste0("~/downloads/", id))
      paired = str_detect(id, "_2")
      rl = as.numeric(str_replace(id, "_.*", ""))
      tab = read.table("sim_info.txt", header=TRUE)
      tab$t_id = str_extract(tab$transcript_id, "ENST.*")
      load("kallisto_reads.rda")
      kallisto_sub = counts_kallisto[match(tab$t_id, rownames(counts_kallisto)),]
      load(paste0("reads_weight1-", rl, ".rda"))
      reads = assays(rse)$counts
      vars = assays(rse)$se^2
      us_sub = apply(reads[match(tab$t_id, rownames(reads)),], 2, as.numeric)
      us_sub[is.na(us_sub)] = 0
      var_sub = apply(vars[match(tab$t_id, rownames(vars)),], 2, as.numeric)
      var_sub[is.na(var_sub)] = Inf
      se_sub = sqrt(var_sub)
      load("hisat2_cov-e.rda")
      st_sub = info[match(tab$t_id, info[,1]),2, drop=FALSE]
      st_sub = apply(st_sub, 2, as.numeric); wids = as.numeric(as.character(tx_wids[match(tab$t_id, tx_wids[,1]), 2])); st_sub = st_sub*wids/rl
      load("rsem_reads.rda")
      rsem_sub = info[match(tab$t_id, rownames(info)),,drop=F]
      load("salmon_reads.rda")
      salmon_sub = salmon[match(tab$t_id, salmon[,1]), 2, drop=F]
      load("hisat_cufflinks_cov.rda")
      cufflinks_sub = info[match(tab$t_id, info[,1]), 2, drop=F]
      cufflinks_sub = apply(cufflinks_sub, 2, as.numeric); wids = as.numeric(as.character(tx_wids[match(tab$t_id, tx_wids[,1]), 2])); cufflinks_sub = cufflinks_sub*wids/rl

      if(paired==T){
            us_sub = us_sub/2
            st_sub = st_sub/2
            cufflinks_sub = cufflinks_sub/2
      }
      real = round(tab$reads)
      bord = 12

      us_err = us_sub[,1] - real
      kallisto_err = kallisto_sub[,1] - real
      st_err = st_sub[,1] - real
      rsem_err = rsem_sub[,1] - real
      salmon_err = salmon_sub[,1] - real
      cufflinks_err = cufflinks_sub[,1]-real

      tss = sum((real-mean(real))^2)
      rss_us = sum((us_err-mean(us_err))^2)
      rss_kallisto = sum((kallisto_err-mean(kallisto_err))^2)
      rss_st = sum((st_err-mean(st_err))^2)
      rss_rsem = sum((rsem_err-mean(rsem_err))^2)
      rss_salmon = sum((salmon_err-mean(salmon_err))^2)
      rss_cufflinks = sum((cufflinks_err-mean(cufflinks_err, na.rm=T))^2, na.rm=T)

      r2_us = 1-rss_us/tss
      r2_kallisto = 1-rss_kallisto/tss
      r2_st = 1-rss_st/tss
      r2_rsem = 1-rss_rsem/tss
      r2_salmon = 1-rss_salmon/tss
      r2_cufflinks = 1-rss_cufflinks/tss

      rmse_us = signif(sqrt(mean(us_err^2, na.rm=T)), 3)
      rmse_kallisto = signif(sqrt(mean(kallisto_err^2, na.rm=T)), 3)
      rmse_st = signif(sqrt(mean(st_err^2, na.rm=T)), 3)
      rmse_rsem = signif(sqrt(mean(rsem_err^2, na.rm=T)), 3)
      rmse_salmon = signif(sqrt(mean(salmon_err^2, na.rm=T)), 3)
      rmse_cufflinks= signif(sqrt(mean(cufflinks_err^2, na.rm=T)), 3)

      png(paste0(id, "_paper.png"), width=600, height=800, units = "px", pointsize=18)
      layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 2:25), nrow=8, byrow=F), heights=c(1, 4, 4, 4, 4, 4, 4, 1), widths=c(1.5, 4, 4, 1.5))
      par(mar=rep(0,4))
      plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
      end = "Single End"
      if(paired==T){end="Paired End"}
      text(x=0.4, y=0.5, srt=90, paste0(rl, "bp ", end, " Simulation Results"), cex=2, pos=3)
      plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
      text(x=0.5, y=0.8, pos=1, "Tx Level", cex=2)

      A = log(us_sub[,1]+1, 2) + log(real+1, 2)
      M = log(us_sub[,1]+1, 2) - log(real+1, 2)
     
      plot(M~A, ylim=c(-bord, bord), xlim=c(0.5, bord)*2, col=rgb(0,0,0,0.1), main = paste0("recountNNLS \nRMSE: ", rmse_us), xlab="MA Plot: Truth - LM",
      xaxt="n", yaxt="n", pch=19, cex=cs); abline(h=0, col=2)
      abline(0, 1, lty=2, col=2); abline(0, -1, lty=2, col=2)
      axis(side=2, at = log(c(1/1000, 1/100, 1/10, 1, 10, 100, 1000),2), labels=parse(text = c('10^-3', '10^-2', '10^-1', '1', '10', '10^2', '10^3')) , las=2)

      A = log(kallisto_sub[,1]+1, 2) + log(real+1, 2)
      M = log(kallisto_sub[,1]+1, 2) - log(real+1, 2)
      plot(M~A, ylim=c(-bord, bord), xlim=c(0.5, bord)*2, col=rgb(0,0,0,0.05), main = paste0("Kallisto \nRMSE: ", rmse_kallisto), xlab="MA Plot Truth - Kallisto",
           xaxt="n", yaxt="n", pch=19, cex=cs); abline(h=0, col=2)
      abline(0, 1, lty=2, col=2); abline(0, -1, lty=2, col=2)
      axis(side=2, at = log(c(1/1000, 1/100, 1/10, 1, 10, 100, 1000),2), labels=parse(text = c('10^-3', '10^-2', '10^-1', '1', '10', '10^2', '10^3')) , las=2)

      A = log(st_sub[,1]+1, 2) + log(real+1, 2)
      M = log(st_sub[,1]+1, 2) - log(real+1, 2)
      plot(M~A, ylim=c(-bord, bord), xlim=c(0.5, bord)*2,  col=rgb(0,0,0,0.05),main = paste0("HISAT2-StringTie \nRMSE: ", rmse_st), xlab="MA Plot Truth - Hisat-Stringtie",
           xaxt="n", yaxt="n", pch=19, cex=cs); abline(h=0, col=2)
      abline(0, 1, lty=2, col=2); abline(0, -1, lty=2, col=2)
      axis(side=2, at = log(c(1/1000, 1/100, 1/10, 1, 10, 100, 1000),2), labels=parse(text = c('10^-3', '10^-2', '10^-1', '1', '10', '10^2', '10^3')) , las=2)

      A = log(cufflinks_sub[,1]+1, 2) + log(real+1, 2)
      M = log(cufflinks_sub[,1]+1, 2) - log(real+1, 2)
      plot(M~A, ylim=c(-bord, bord), xlim=c(0.5, bord)*2,  col=rgb(0,0,0,0.05),main = paste0("HISAT2-Cufflinks \nRMSE: ", rmse_cufflinks), xlab="MA Plot Truth - cufflinks",
           xaxt="n", yaxt="n", pch=19, cex=cs); abline(h=0, col=2)
      abline(0, 1, lty=2, col=2); abline(0, -1, lty=2, col=2)
      axis(side=2, at = log(c(1/1000, 1/100, 1/10, 1, 10, 100, 1000),2), labels=parse(text = c('10^-3', '10^-2', '10^-1', '1', '10', '10^2', '10^3')) , las=2)

      A = log(rsem_sub[,1]+1, 2) + log(real+1, 2)
      M = log(rsem_sub[,1]+1, 2) - log(real+1, 2)
      plot(M~A, ylim=c(-bord, bord), xlim=c(0.5, bord)*2,  col=rgb(0,0,0,0.05),main = paste0("RSEM \nRMSE: ", rmse_rsem), xlab="MA Plot Truth - RSEM",
           xaxt="n", yaxt="n", pch=19, cex=cs); abline(h=0, col=2)
      abline(0, 1, lty=2, col=2); abline(0, -1, lty=2, col=2)
      axis(side=2, at = log(c(1/1000, 1/100, 1/10, 1, 10, 100, 1000),2), labels=parse(text = c('10^-3', '10^-2', '10^-1', '1', '10', '10^2', '10^3')) , las=2)

      A = log(salmon_sub[,1]+1, 2) + log(real+1, 2)
      M = log(salmon_sub[,1]+1, 2) - log(real+1, 2)
      plot(M~A, ylim=c(-bord, bord), xlim=c(0.5, bord)*2,  col=rgb(0,0,0,0.05),main = paste0("Salmon \nRMSE: ", rmse_salmon), xlab="MA Plot Truth - salmon",
           xaxt="n", yaxt="n", pch=19, cex=cs); abline(h=0, col=2)
      abline(0, 1, lty=2, col=2); abline(0, -1, lty=2, col=2)
      axis(side=2, at = log(c(1/1000, 1/100, 1/10, 1, 10, 100, 1000),2), labels=parse(text = c('10^-3', '10^-2', '10^-1', '1', '10', '10^2', '10^3')) , las=2)
      axis(side=1, at = log(c(1, 10, 100, 1000, 10000, 100000, 1000000), 2), 
            labels=parse(text = c(1, 10, 100, '10^3', '10^4', '10^5', '10^6')))

      plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
      plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
      text(x=0.5, y=0.8, pos=1, "Gene Level", cex=2)
      
      us_sub = cbind(by(us_sub[,1], tab$gene_id, sum))
      kallisto_sub = cbind(by(kallisto_sub[,1], tab$gene_id, sum))
      st_sub = cbind(by(st_sub[,1], tab$gene_id, sum))
      rsem_sub = cbind(by(rsem_sub[,1], tab$gene_id, sum))
      salmon_sub = cbind(by(salmon_sub[,1], tab$gene_id, sum))
      cufflinks_sub = cbind(by(cufflinks_sub[,1], tab$gene_id, sum))
      real = by(real, tab$gene_id, sum)
      bord = 16.5

      us_err = us_sub[,1] - real
      kallisto_err = kallisto_sub[,1] - real
      st_err = st_sub[,1] - real
      rsem_err = rsem_sub[,1] - real
      salmon_err = salmon_sub[,1] - real
      cufflinks_err = cufflinks_sub[,1] - real

      rmse_us = signif(sqrt(mean(us_err^2, na.rm=T)), 3)
      rmse_kallisto = signif(sqrt(mean(kallisto_err^2, na.rm=T)), 3)
      rmse_st = signif(sqrt(mean(st_err^2, na.rm=T)), 3)
      rmse_rsem = signif(sqrt(mean(rsem_err^2, na.rm=T)), 3)
      rmse_salmon = signif(sqrt(mean(salmon_err^2, na.rm=T)), 3)
      rmse_cufflinks = signif(sqrt(mean(cufflinks_err^2, na.rm=T)), 3)

      tss = sum((real-mean(real))^2)
      rss_us = sum((us_err-mean(us_err))^2)
      rss_kallisto = sum((kallisto_err-mean(kallisto_err))^2)
      rss_st = sum((st_err-mean(st_err))^2)
      rss_rsem = sum((rsem_err-mean(rsem_err))^2)
      rss_salmon = sum((salmon_err-mean(salmon_err))^2)
      rss_cufflinks = sum((cufflinks_err-mean(cufflinks_err, na.rm=T))^2, na.rm=T)

      r2_us = 1-rss_us/tss
      r2_kallisto = 1-rss_kallisto/tss
      r2_st = 1-rss_st/tss
      r2_rsem = 1-rss_rsem/tss
      r2_salmon = 1-rss_salmon/tss
      r2_cufflinks = 1-rss_cufflinks/tss

      A = log(us_sub[,1]+1, 2) + log(real+1, 2)
      M = log(us_sub[,1]+1, 2) - log(real+1, 2)
      plot(M~A, ylim=c(-bord, bord), xlim=c(0.5, bord)*2, col=rgb(0,0,0,0.1), main = paste0("recountNNLS \nRMSE: ", rmse_us), xlab="MA Plot: Truth - LM",
           xaxt="n", yaxt="n", pch=19, cex=cs); abline(h=0, col=2)
      abline(0, 1, lty=2, col=2); abline(0, -1, lty=2, col=2)
      axis(side=4, at = log(c(1/10000, 1/1000, 1/100, 1/10, 1, 10, 100, 1000, 10000),2), 
            labels=parse(text = c('10^-4', '10^-3', '10^-2', '10^-1', '1', '10', '10^2', '10^3', '10^4')) , las=2)

      A = log(kallisto_sub[,1]+1, 2) + log(real+1, 2)
      M = log(kallisto_sub[,1]+1, 2) - log(real+1, 2)
      plot(M~A, ylim=c(-bord, bord), xlim=c(0.5, bord)*2, col=rgb(0,0,0,0.1), main = paste0("Kallisto \nRMSE: ", rmse_kallisto), xlab="MA Plot Truth - Kallisto",
           xaxt="n", yaxt="n", pch=19, cex=cs); abline(h=0, col=2)
      abline(0, 1, lty=2, col=2); abline(0, -1, lty=2, col=2)
      axis(side=4, at = log(c(1/10000, 1/1000, 1/100, 1/10, 1, 10, 100, 1000, 10000),2), 
            labels=parse(text = c('10^-4', '10^-3', '10^-2', '10^-1', '1', '10', '10^2', '10^3', '10^4')) , las=2)

      A = log(st_sub[,1]+1, 2) + log(real+1, 2)
      M = log(st_sub[,1]+1, 2) - log(real+1, 2)
      plot(M~A, ylim=c(-bord, bord), xlim=c(0.5, bord)*2,  col=rgb(0,0,0,0.1),main = paste0("HISAT2-StringTie \nRMSE: ", rmse_st), xlab="MA Plot Truth - Hisat-Stringtie",
           xaxt="n", yaxt="n", pch=19, cex=cs); abline(h=0, col=2)
      abline(0, 1, lty=2, col=2); abline(0, -1, lty=2, col=2)
      axis(side=4, at = log(c(1/10000, 1/1000, 1/100, 1/10, 1, 10, 100, 1000, 10000),2), 
            labels=parse(text = c('10^-4', '10^-3', '10^-2', '10^-1', '1', '10', '10^2', '10^3', '10^4')) , las=2)

      A = log(cufflinks_sub[,1]+1, 2) + log(real+1, 2)
      M = log(cufflinks_sub[,1]+1, 2) - log(real+1, 2)
      plot(M~A, ylim=c(-bord, bord), xlim=c(0.5, bord)*2,  col=rgb(0,0,0,0.1),main = paste0("HISAT2-Cufflinks \nRMSE: ", rmse_cufflinks), xlab="MA Plot Truth - RSEM",
           xaxt="n", yaxt="n", pch=19, cex=cs); abline(h=0, col=2)
      abline(0, 1, lty=2, col=2); abline(0, -1, lty=2, col=2)
      axis(side=4, at = log(c(1/10000, 1/1000, 1/100, 1/10, 1, 10, 100, 1000, 10000),2), 
            labels=parse(text = c('10^-4', '10^-3', '10^-2', '10^-1', '1', '10', '10^2', '10^3', '10^4')) , las=2)

      A = log(rsem_sub[,1]+1, 2) + log(real+1, 2)
      M = log(rsem_sub[,1]+1, 2) - log(real+1, 2)
      plot(M~A, ylim=c(-bord, bord), xlim=c(0.5, bord)*2,  col=rgb(0,0,0,0.1),main = paste0("RSEM \nRMSE: ", rmse_rsem), xlab="MA Plot Truth - RSEM",
           xaxt="n", yaxt="n", pch=19, cex=cs); abline(h=0, col=2)
      abline(0, 1, lty=2, col=2); abline(0, -1, lty=2, col=2)
      axis(side=4, at = log(c(1/10000, 1/1000, 1/100, 1/10, 1, 10, 100, 1000, 10000),2), 
            labels=parse(text = c('10^-4', '10^-3', '10^-2', '10^-1', '1', '10', '10^2', '10^3', '10^4')) , las=2)

      A = log(salmon_sub[,1]+1, 2) + log(real+1, 2)
      M = log(salmon_sub[,1]+1, 2) - log(real+1, 2)
      plot(M~A, ylim=c(-bord, bord), xlim=c(0.5, bord)*2,  col=rgb(0,0,0,0.1),main = paste0("Salmon \nRMSE: ", rmse_salmon), xlab="MA Plot Truth - RSEM",
           xaxt="n", yaxt="n", pch=19, cex=cs); abline(h=0, col=2)
      abline(0, 1, lty=2, col=2); abline(0, -1, lty=2, col=2)
      axis(side=4, at = log(c(1/10000, 1/1000, 1/100, 1/10, 1, 10, 100, 1000, 10000),2), 
            labels=parse(text = c('10^-4', '10^-3', '10^-2', '10^-1', '1', '10', '10^2', '10^3', '10^4')) , las=2)
      axis(side=1, at = log(c(10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000), 2), 
            labels=parse(text = c(10, 100, '10^3', '10^4', '10^5', '10^6', '10^7', '10^8', '10^9')))

      par(xpd=NA)
      plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
      plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
      plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
      text(x=0.6, y=0.5, srt=270, "recountNNLS", cex=1.4, pos=3)
      plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
      text(x=0.6, y=0.5, srt=270, "Kallisto", cex=1.4, pos=3)
      plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
      text(x=0.6, y=0.5, srt=270, "H2-StringTie", cex=1.4, pos=3)
      plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
      text(x=0.6, y=0.5, srt=270, "H2-Cufflinks", cex=1.4, pos=3)
      plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
      text(x=0.6, y=0.5, srt=270, "RSEM", cex=1.4, pos=3)
      plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
      text(x=0.6, y=0.5, srt=270, "Salmon", cex=1.4, pos=3)

      par(xpd=F)
      dev.off()
}

