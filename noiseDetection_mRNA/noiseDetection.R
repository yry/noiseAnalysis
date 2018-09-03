###################################
##plot the p2pPCC vs the abundance
###################################

	abn.filename = "matrix.abn"
	pcc.filename = "matrix.pcc"
	pdf.filename = "p2p_correlations.pdf"
#	print(abn.filename)
#	print(pcc.filename)

	abn.matrix <- read.csv(abn.filename,sep=',',header=FALSE)
	pcc.matrix <- read.csv(pcc.filename,sep=',',header=FALSE)
	n.samples=4

	pdf(pdf.filename,width=15,height=5)
	par(mfrow=c(1,3))
	for(i in 1:n.samples)
	{
	   print(i) 
	  abn.bins2 <- floor(2*log2(abn.matrix[,i]+1))
	  boxplot(pcc.matrix[,i]~abn.bins2,main=paste("sample",i), 
		  ylim=c(-1,1),xlab='Abundance', ylab='p2p PCC',xaxt='n')
	  abline(v=10,col='red'); abline(v=12,col='blue')
	  abline(h=0.7,col='red');abline(h=0.6,col='blue')
	  xlab <- (1:40)/2
	  axis(1,at=1:40,labels=xlab,las=2)
	}
	dev.off()

###################################
##determine the singa/noise thresholds
###################################

	abn.filename = "matrix.abn"
	pcc.filename = "matrix.pcc"
	csv.filename = "noiseThresholds.csv"
	print(abn.filename)
	print(pcc.filename)
	print(csv.filename)

	abn.matrix <- read.csv(abn.filename,sep=',',header=FALSE)
	pcc.matrix <- read.csv(pcc.filename,sep=',',header=FALSE)
	n.samples=4
	
	iqr.thr = 0.6
	noise.thrU <- rep(-1,n.samples) ##upper noise bound
	noise.thrL <- rep(-1,n.samples) ##lower noise bound
	for(i in 1:n.samples)
	{
	  abn.bins <- floor(log2(abn.matrix[,i]+1))
	  pcc.summary <- tapply(pcc.matrix[,i],as.factor(abn.bins),FUN=summary)
	  for(j in 2:length(pcc.summary))
	  {
		##median above the thr
		##starting from the left (min abundance)
		if(pcc.summary[[j]][3] > iqr.thr)
		{
		noise.thrL[i] = j
		break
		}
	  }
	  for(j in seq(length(pcc.summary)-1,2,-1))
	  {
		##median above the thr
		##starting from the right (max abundance)
		if(pcc.summary[[j-1]][3] < iqr.thr)
		{
		noise.thrU[i] = j
		break
		}
	  }

	}
	library(MASS)
	noise.thr = cbind(seq(1,n.samples),noise.thrL,noise.thrU)
	colnames(noise.thr)=c("sample","minNoise","maxNoise")
	write.matrix(noise.thr,file=csv.filename,sep=',')


