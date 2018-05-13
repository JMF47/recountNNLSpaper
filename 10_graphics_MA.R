####################################################
### MA Plots (supplements)
####################################################
cs = 0.15
ids = c("37_1", "37_2", "50_1", "50_2", "75_1", "75_2", "100_1", "100_2", "150_1", "150_2")

for(id in ids){
      scenario = str_split_fixed(id, "_", n=2)
      rl = scenario[1]
      paired = scenario[2]
      load(paste0("/dcl01/leek/data/ta_poc/geuvadis/simulation/", id, "/", rl, "_", paired, "_info.rda"))
      ses = recountNNLSse
      scores = recountNNLSscore

      err = info-truth$reads
      mlab = "MAE"
      metric = round(apply(abs(err), 2, mean, na.rm=T), 3)
      
      png(paste0('/dcl01/leek/data/ta_poc/graphics/', id, "_MA.png"), 
            width=400, height=800, units = "px", pointsize=18, type="cairo")
      layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 2:8, 9:15), nrow=7, byrow=F), 
            heights=c(1, 4, 4, 4, 4, 4, 1), widths=c(1, 4, 0.5))

      ### texts
      par(mar=rep(0,4))
      plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
      end = "Single End"
      if(paired==2){end="Paired End"}
      text(x=0.4, y=0.5, srt=90, paste0(rl, "bp ", end, " Simulation Results"), cex=2, pos=3)
      plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
      text(x=0.5, y=0.8, pos=1, "Tx Level", cex=2)

      bord = 12/2

      ### First set of MA plots
      for(i in 1:5){
            A = (log(info[,i]+1, 2) + log(truth$reads+1, 2))/2
            M = log(info[,i]+1, 2) - log(truth$reads+1, 2)

            plot(M~A, ylim=c(-bord*2, bord*2), xlim=c(0.5, bord)*2, col=rgb(0,0,0,0.1), 
                  main = "", xlab="MA Plot: Truth - LM", xaxt="n", yaxt="n", pch=19, cex=cs)
            # abline(h=0, col=2)
            # abline(0, 1, lty=2, col=2); abline(0, -1, lty=2, col=2)
            axis(side=2, at = log(c(1/1000, 1/100, 1/10, 1, 10, 100, 1000),2), 
                  labels=parse(text = c('10^-3', '10^-2', '10^-1', '1', '10', '10^2', '10^3')) , las=2)
            window = par()$usr
            ywid = window[4]-window[3]
            text((window[1]+window[2])/2,  window[3]+ywid/10, labels=paste0(mlab, ": ", metric[i]), cex=1.5)
      }
      axis(side=1, at = log(c(1, 10, 100, 1000, 10000, 100000, 1000000), 2),
            labels=parse(text = c(1, 10, 100, '10^3', '10^4', '10^5', '10^6')))

      par(xpd=NA)
      plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
      plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
      plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
      text(x=0.6, y=0.5, srt=270, "recountNNLS", cex=1.4, pos=3)
      plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
      text(x=0.6, y=0.5, srt=270, "Kallisto", cex=1.4, pos=3)
      plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
      text(x=0.6, y=0.5, srt=270, "H2-Cufflinks", cex=1.4, pos=3)
      plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
      text(x=0.6, y=0.5, srt=270, "RSEM", cex=1.4, pos=3)
      plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
      text(x=0.6, y=0.5, srt=270, "Salmon", cex=1.4, pos=3)

      par(xpd=F)
      dev.off()
}

####################################################
### MA Plots - RSEM-based
####################################################
load("/dcl01/leek/data/ta_poc/geuvadis/simulation/rsem_based/rsem_based_0.rda")
out = cbind(out,out[,4])
err = out-count_mat[,1]
mlab = "MAE"
metric = round(apply(abs(err), 2, mean, na.rm=T), 3)

png(paste0('/dcl01/leek/data/ta_poc/graphics/RSEM_MA.png'), 
      width=400, height=800, units = "px", pointsize=18, type="cairo")
layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 2:8, 9:15), nrow=7, byrow=F), 
      heights=c(1, 4, 4, 4, 4, 4, 1), widths=c(1, 4, 0.5))

### texts
par(mar=rep(0,4))
plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
end = "Single End"
end="Paired End"
text(x=0.4, y=0.5, srt=90, paste0("RSEM-based Simulation Results"), cex=2, pos=3)
plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
text(x=0.5, y=0.8, pos=1, "Tx Level", cex=2)

bord = 10
cs = 0.4
for(i in 1:5){
      color = rgb(0,0,0,0.1)
      if(i==4) color = rgb(0, 0, 0, 0)
      A = (log(out[,i]+1, 2) + log(count_mat[,1]+1, 2))/2
      M = log(out[,i]+1, 2) - log(count_mat[,1]+1, 2)

      plot(M~A, ylim=c(-10, 10), xlim=c(0, 18), col=color, 
            main = "", xlab="MA Plot: Truth - LM", xaxt="n", yaxt="n", pch=19, cex=cs)
      axis(side=2, at = log(c(1/1000, 1/100, 1/10, 1, 10, 100, 1000),2), 
            labels=parse(text = c('10^-3', '10^-2', '10^-1', '1', '10', '10^2', '10^3')) , las=2)
      window = par()$usr
      ywid = window[4]-window[3]
      if(i!=4){text((window[1]+window[2])/2,  window[3]+ywid/10, labels=paste0(mlab, ": ", metric[i]), cex=1.5)}
      
}
axis(side=1, at = log(c(1, 100, 10000, 10^8, 10^16), 2),
      labels=parse(text = c(1, 100, '10^4', '10^8', '10^16')))

par(xpd=NA)
plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
text(x=0.6, y=0.5, srt=270, "recountNNLS", cex=1.4, pos=3)
plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
text(x=0.6, y=0.5, srt=270, "Kallisto", cex=1.4, pos=3)
plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
text(x=0.6, y=0.5, srt=270, "H2-Cufflinks", cex=1.4, pos=3)
plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
text(x=0.6, y=0.5, srt=270, "RSEM", cex=1.4, pos=3)
plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
text(x=0.6, y=0.5, srt=270, "Salmon", cex=1.4, pos=3)

par(xpd=F)
dev.off()