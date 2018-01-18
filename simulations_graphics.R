	##########################################################################################
	### Look at errors
	##########################################################################################
	err = info - truth$transcriptsReads
	err2 = err^2
		err = l2info - l2truth
		err2 = err^2
	sqrt(apply(err2, 2, mean))

	err_gene = info_gene - truth_gene
	err2_gene = err_gene^2
		err_gene = l2info_gene - l2truth_gene
		err2_gene = err_gene^2
	sqrt(apply(err2_gene, 2, mean))

	### Look at number of called expressed
	apply(l2info>0, 2, sum)
	sum(l2truth>0)

	apply(l2info_gene>0, 2, sum)
	sum(l2truth_gene>0)
	
	### Correlations
	cor(info)
	cor(l2info)

	apply(l2info, 2, cor, l2truth, method="spearman")
	apply(l2info_gene, 2, cor, l2truth_gene, method="spearman")

	##########################################################################################
	### Custom MA plots
	##########################################################################################
	M = l2info-l2truth; A = l2info+l2truth; M = round(M*4)/4; A = round(A*4)/4
	xs = seq(0, 20, by=0.25); ys = seq(-10, 10, by=0.25)

	cells = NULL;coordinates = NULL
	for(i in 1:6){
		tmp = NULL
		for(j in 1:length(xs)){
			for(k in 1:length(ys)){
				message(i, "-", j, "-", k)	
				tmp = c(tmp, sum(A[,i] == xs[j] & M[,i]==ys[k] ))
				coordinates = rbind(coordinates, c(xs[j], ys[k]))
			}
		}
		cells = cbind(cells, tmp)
	}
	zz = cells[41,] ### Setting (0,0) to have 0 counts so it doesnt skew the colors
	cells[41,] = rep(0, 6)
	lcells = log(cells+1, 2)

	colfunc <- colorRampPalette(c("white", "black"))
	cols = colfunc(101)

	color_scale = round(lcells/max(lcells)*100)
	
	layout(matrix(c(1:9), nrow=9, byrow=F), heights=c(1, 4, 4, 4, 4, 4, 4, 1, 2), widths=1)
	par(mar=c(0, 0, 0, 0))
	plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
    text(x=0.5, y=0.8, pos=1, "Tx Level", cex=2)

	for(i in 1:6){
		par(mar=c(0, 2, 0, 0))
		plot(c(min(xs), max(xs)), c(min(ys), max(ys)), ty="n", xlab="", ylab="", xaxt="n", yaxt="n")
		axis(side=2, at=c(-8, -4, 0, 4, 8))
		abline(h=0, col=2)
		rect(xleft=coordinates[,1], xright=coordinates[,1]+0.25, ybot = coordinates[,2], ytop = coordinates[,2]+0.25,
			col = cols[color_scale[,i]+1], border=cols[color_scale[,i]+1])	
	}
	plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
	plot(c(-10, 110), c(0, 1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
	par(mar=c(3,0,0,0))
	rect(xleft=0:100, xright=1:101, ytop=0.8, ybot=0.2, col=cols, border=cols)
		positions = c(0, 25, 50, 75, 100); p = positions[-1]/100; counts = round(2^(max(log(cells+1,2))*p))
		axis(side=1, at=positions, labels=parse(text = c(0, counts)))




	ids = paste0(M, "-", A)




