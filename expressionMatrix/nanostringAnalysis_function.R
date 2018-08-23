##set of functions for the nanostring analysis
##load this file at the begining of the analysis

#' Function to calcualte the correlation matrices using a sliding window
#'
#' This function generates an average correlation coefficient for every sliding window
#' for each sample, in the expression matrix. 
#' @param expression.matrix expression matrix, normalized
#' @param method one of the correlation coefficients, Pearson, Spearman or Kendall
#' @param n.elements.per.window number of elements to have in a window, default 40
#' @keywords correlation coefficients on sliding windows
#' @export
#' @examples
#' calculateCorrelationMatrices(expression.matrix = expr.mat.host, method="pearson")
#' calculateCorrelationMatrices(expression.matrix = expr.mat.host, method="spearman")
#' calculateCorrelationMatrices(expression.matrix = expr.mat.host, method="kendall")
calculateCorrelationMatrices <- function(expression.matrix,
                                         method = c("pearson", "kendall", "spearman"),n.elements.per.window=40)
{
  cat("the input matrix has ",nrow(expression.matrix)," rows and ",ncol(expression.matrix)," cols\n")
  cat("number of genes: ",nrow(expression.matrix),"\n")
  cat("number of samples: ",ncol(expression.matrix),"\n")
  
  if(!exists("n.elements.per.window"))
  {
    n.elements.per.window = floor(nrow(expression.matrix) / 10)
    print("calculating the number of elements per window")
  }
  cat("the number of elements per window is ",n.elements.per.window,"\n")
  
  ##the matrix contains the correlation coefficient (CC) for each sample, at each abundance level.
  ##the abundance level is maintained in abn.matrix variable
  corr.matrix = rep(0,ncol(expression.matrix)*(nrow(expression.matrix)-n.elements.per.window+1))
  corr.matrix = matrix(corr.matrix,ncol=ncol(expression.matrix))
  abn.matrix = rep(0,ncol(expression.matrix)*(nrow(expression.matrix)-n.elements.per.window+1))
  abn.matrix = matrix(abn.matrix,ncol=ncol(expression.matrix))
  
  for(j in 1:ncol(expression.matrix))
  {
    cat("working with sample ",j,"\n")
    sorted.matrix = expression.matrix[order(expression.matrix[,j]),]
    for(idx in 1:(nrow(expression.matrix)-n.elements.per.window+1))
    {
      #focus on a sliding window, initialize the correlation vector
      correlation.vector = rep(0,ncol(expression.matrix))
      for(k in 1:ncol(expression.matrix))
      {
        if(j != k)
        {
          col.j = sorted.matrix[idx:(idx+n.elements.per.window-1),j];names(col.j)=rep("",length(col.j))
          col.k = sorted.matrix[idx:(idx+n.elements.per.window-1),k];names(col.k)=rep("",length(col.k))
          cor.in = cbind(col.j,col.k)
          if(method == "pearson"){correlation.vector[k] = cor(cor.in,method="pearson")[1,2]}
          if(method == "spearman"){correlation.vector[k] = cor(cor.in,method="spearman")[1,2]}
          if(method == "kendall"){correlation.vector[k] = cor(cor.in,method="kendall")[1,2]}
        }#end if
      }#end for k
      corr.matrix[idx,j] = mean(correlation.vector[correlation.vector != 0])
      abn.matrix[idx,j] = mean(sorted.matrix[idx:(idx+n.elements.per.window-1),j])
    }#end for on the sliding windows
  }
  
  returnObject <- list("abn" = abn.matrix, "corr" = corr.matrix)
  return(returnObject)
}

##optimize the number of elemens per window
##the optimal size window is the one for which no changes are observed.
##specify the column/sample for which the sarch is performed
##the output is the recommended minimum number of elements per window
optimize.elements.per.window <- function(expression.matrix, selected.col = 1,
                                         method = c("pearson", "kendall", "spearman"))
{
  n.elements.per.window = floor(nrow(expression.matrix) / 10)
  n.elts = 10
  k = 1; #projection of the search space
  #selected correlations and the selection of corresponding abundance
  P.corr.out <- matrix(rep(-1,nrow(expression.matrix)*50),ncol=50)
  P.abn.out  <- matrix(rep(-1,nrow(expression.matrix)*50),ncol=50)
  
  while(n.elts < 2 * n.elements.per.window)
  {
    print(n.elts)
    P.corr.k <- calculateCorrelationMatrices(expression.matrix = expr.mat.host, 
                                             method=method,n.elements.per.window = n.elts)
    for(i in 1:nrow(P.corr.k$corr))
    {
      P.corr.out[i,k] = P.corr.k$corr[i,selected.col]
      P.abn.out[i,k]  = P.corr.k$abn[i,selected.col]
    }
    #plot.PSK(paste("P",k,".pdf",sep=""),abn.matrix = P.corr.k$abn,P.matrix = P.corr.k$corr)
    k = k + 1
    n.elts = n.elts + 5
  }
  ##plot the entries != -1 [initialization value]
  ##determine the stability value
  ##the plotting is done relative to the average abundance in the window
  pdf("paramSearch.pdf",width = 20,height=8)
  par(mfrow=c(2,4))
  for(k.idx in (1:(k-1)))
  {
    plot.y <- P.corr.out[P.corr.out[,k.idx]!=-1,k.idx]
    plot.x <- P.abn.out[P.abn.out[,k.idx]!=-1,k.idx]
    winLen = 10+(k.idx-1)*5
    plot(plot.x,plot.y,#col=k.idx,
         type='l',ylim=c(-1,1),xlim=c(0,max(P.abn.out)),main=paste("win len:",winLen));
  }
  par(new=FALSE)
  dev.off()
  
  ##determine the minimum number of elements per window.
  ##compare the stability of the output plots
  for(i in (2:(k-1)))
  {
    maxN <- sum(P.corr.out[,i] != -1)
    corr.in  <- P.corr.out[1:maxN,c(i-1,i)]
    corr.in[is.na(corr.in)] = 0
    corr.win <- cor(corr.in, method = "pearson")[1,2]
    winLen = 10+(i-1)*5
    print(paste(winLen,corr.win))
    if(corr.win > 0.95)
    {
      return(winLen)
    }
  }
  
}

##this function removes the noise reads from expression matrix
##uses as input a correlation matrix and teh corresponding abundance matrix
##all entries below the noise threshold are replaced with the noise threshold
## corr.thr is the correlation threshold to be applied
removeNoise.correlation <- function(expression.matrix, CC.matrix, abn.matrix, corr.thr)
{
  ##the 4 rows for these matrices correspond to the non-smoothed and smoothed analysis (10,25,50)
  CC.thr  = matrix(rep(0,4*ncol(expression.matrix)),nrow=4)
  abn.thr = matrix(rep(0,4*ncol(expression.matrix)),nrow=4)
  
  ##initialize the matrix which will be returned with the original matrix
  expression.matrix.noNoise = expression.matrix
  print("Sample,AbnThrRaw, AbnThrLoess10,AbnThrLoess25,AbnThrLoess50")
  for(j in 1:ncol(expression.matrix))
  {
    loessMod10 <- loess(CC.matrix[,j] ~ abn.matrix[,j], span=0.10)
    smoothed10 <- predict(loessMod10)
    loessMod25 <- loess(CC.matrix[,j] ~ abn.matrix[,j], span=0.25)
    smoothed25 <- predict(loessMod25)
    loessMod50 <- loess(CC.matrix[,j] ~ abn.matrix[,j], span=0.50)
    smoothed50 <- predict(loessMod50)
    
    CC.thr.raw = CC.matrix[,j] > corr.thr
    CC.thr.10  = smoothed10 > corr.thr
    CC.thr.25  = smoothed25 > corr.thr
    CC.thr.50  = smoothed50 > corr.thr
    
    ##the indexes for which the transition is achieved
    CC.thr[1,j] = max(which(CC.thr.raw == FALSE))
    CC.thr[2,j] = max(which(CC.thr.10  == FALSE))
    CC.thr[3,j] = max(which(CC.thr.25  == FALSE))
    CC.thr[4,j] = max(which(CC.thr.50 == FALSE))
    
    ##the correspondong abundances
    abn.thr[1,j] = abn.matrix[CC.thr[1,j],j]
    abn.thr[2,j] = abn.matrix[CC.thr[2,j],j]
    abn.thr[3,j] = abn.matrix[CC.thr[3,j],j]
    abn.thr[4,j] = abn.matrix[CC.thr[4,j],j]
    
    expression.matrix.noNoise[expression.matrix[,j]<abn.thr[1,j],j] = abn.thr[1,j]
    abn.thr[1,j] = floor(abn.thr[1,j])
    abn.thr[2,j] = floor(abn.thr[2,j])
    abn.thr[3,j] = floor(abn.thr[3,j])
    abn.thr[4,j] = floor(abn.thr[4,j])
    print(paste(j,abn.thr[1,j],abn.thr[2,j],abn.thr[3,j],abn.thr[4,j],sep=','))
  }
  # write.csv(expression.matrix.noNoise,
  # file='onlyHost_notNorm_expressionMatrix_noNoise.csv')
  return(expression.matrix.noNoise)
}

#' Function to calculate the Kullback-Leibler divergence matrices using a sliding window
#'
#' This function generates an average KLdivergence value for every sliding window
#' for each sample, in the expression matrix. 
#' @param expression.matrix expression matrix, normalized
#' @keywords KLdivergences on sliding windows
#' @export
#' @examples
calculateKLMatrix <- function(expression.matrix)
{
  cat("the input matrix has ",nrow(expression.matrix)," rows and ",ncol(expression.matrix)," cols\n")
  cat("number of genes: ",nrow(expression.matrix),"\n")
  cat("number of samples: ",ncol(expression.matrix),"\n")
  
  n.elements.per.window = floor(nrow(expression.matrix) / 10)
  cat("the number of elements per window is ",n.elements.per.window,"\n")
  n.classes.KL = 10 ##number of classes for the KL divergence
  
  KL.matrix = rep(0,ncol(expression.matrix)*(nrow(expression.matrix)-n.elements.per.window+1))
  KL.matrix = matrix(KL.matrix,ncol=ncol(expression.matrix))
  abn.matrix = rep(0,ncol(expression.matrix)*(nrow(expression.matrix)-n.elements.per.window+1))
  abn.matrix = matrix(abn.matrix,ncol=ncol(expression.matrix))
  
  #calculate the KL div on the standardized values
  for(j in 1:ncol(expression.matrix))
  {
    cat("working with sample ",j,"\n")
    sorted.matrix = expression.matrix[order(expression.matrix[,j]),]
    for(idx in 1:(nrow(expression.matrix)-n.elements.per.window+1))
    {
      ##calculate the p distribution (corresponding to col j) once
      col.j = sorted.matrix[idx:(idx+n.elements.per.window-1),j]
      mean.j = mean(col.j); sd.j = sd(col.j)
      z.j = (col.j - mean.j)/sd.j; names(z.j) = rep("",length(z.j))
      min.j = min(z.j); max.j = max(z.j); 
      
      #focus on a sliding window
      KL.div.vector = rep(0,ncol(expression.matrix))
      for(k in 1:ncol(expression.matrix))
      {
        if(j != k)
        {
          col.k = sorted.matrix[idx:(idx+n.elements.per.window-1),k]
          mean.k = mean(col.k); sd.k = sd(col.k)
          z.k = (col.k - mean.k)/sd.k; names(z.k) = rep("",length(z.k))
          min.k = min(z.k); max.k = max(z.k); 
          
          min.all = ifelse(min.j<min.k,min.j,min.k)
          max.all = ifelse(max.j>max.k,max.j,max.k)
          step.KL = (max.all - min.all)/n.classes.KL
          
          ##calculate the p distribution
          p.dist = rep(0,n.classes.KL+1)
          for(idx1 in 1:length(z.j))
          {
            KL.class = floor((z.j[idx1]-min.all)/step.KL) + 1
            p.dist[KL.class] = p.dist[KL.class]+1
          }#end for for the p distirbution
          p.dist = p.dist/n.elements.per.window
          
          ##calculate the q distribution
          q.dist = rep(0,n.classes.KL+1)
          for(idx2 in 1:length(z.k))
          {
            KL.class = floor((z.k[idx2]-min.all)/step.KL) + 1
            q.dist[KL.class] = q.dist[KL.class]+1
          }#end for for the q distribution
          q.dist = q.dist/n.elements.per.window
          
          #to avoid invalid fractions and logs, the comparison of distirbutions
          #can be done on the 1 - dist
          p.dist = 1 - p.dist; q.dist = 1 - q.dist
          KL.calculate = p.dist*log2(p.dist/q.dist)
          KL.div = sum(KL.calculate)
          KL.div.vector[k] = KL.div
        }#end if j != k
      }
      ##calculate the average KL div across samples
      ##maintain a variable with the average abundance in the selected window
      KL.matrix[idx,j] = mean(KL.div.vector[KL.div.vector != 0])
      abn.matrix[idx,j] = mean(sorted.matrix[idx:(idx+n.elements.per.window-1),j])
    }
  }
  
  returnObject <- list("abn" = abn.matrix, "KL.matrix" = KL.matrix)
  return(returnObject)
}

#' Function to calculate the Kolmogorov-Smirnov matrices using a sliding window
#' d	   =    	max[abs{S1(Y)-S2(Y)}]
#'
#' This function generates an average KS value for every sliding window
#' for each sample, in the expression matrix. 
#' @param expression.matrix expression matrix, normalized
#' @keywords KS values on sliding windows
#' @export
#' @examples
calculateKSMatrix <- function(expression.matrix)
{
  cat("the input matrix has ",nrow(expression.matrix)," rows and ",ncol(expression.matrix)," cols\n")
  cat("number of genes: ",nrow(expression.matrix),"\n")
  cat("number of samples: ",ncol(expression.matrix),"\n")
  
  n.elements.per.window = floor(nrow(expression.matrix) / 10)
  cat("the number of elements per window is ",n.elements.per.window,"\n")
  
  KS.matrix = rep(0,ncol(expression.matrix)*(nrow(expression.matrix)-n.elements.per.window+1))
  KS.matrix = matrix(KS.matrix,ncol=ncol(expression.matrix))
  abn.matrix = rep(0,ncol(expression.matrix)*(nrow(expression.matrix)-n.elements.per.window+1))
  abn.matrix = matrix(abn.matrix,ncol=ncol(expression.matrix))
  
  #calculate the KL div on the standardized values
  for(j in 1:ncol(expression.matrix))
  {
    cat("working with sample ",j,"\n")
    sorted.matrix = expression.matrix[order(expression.matrix[,j]),]
    for(idx in 1:(nrow(expression.matrix)-n.elements.per.window+1))
    {
      col.j = sorted.matrix[idx:(idx+n.elements.per.window-1),j]
      names(col.j) = rep("",length(col.j))
      
      #focus on a sliding window
      KS.vector = rep(0,ncol(expression.matrix))
      for(k in 1:ncol(expression.matrix))
      {
        if(j != k)
        {
          col.k = sorted.matrix[idx:(idx+n.elements.per.window-1),k]
          names(col.k) = rep("",length(col.k))
          
          col.j.sorted = col.j[order(col.j)]
          col.k.sorted = col.k[order(col.k)]
          ##transform these values to cumulative distributions
          ##then apply the KS statistic
          col.j.cumulative = col.j.sorted/sum(col.j.sorted)
          col.k.cumulative = col.k.sorted/sum(col.k.sorted)
          KS.dist = rep(0,length(col.j))
          for(l in 1:length(col.j))
          {
            KS.dist[l]=abs(sum(col.j.cumulative[1:l]) - sum(col.k.cumulative[1:l]))
          }
          KS.vector[k] = max(KS.dist)
        }#end if j != k
      }
      ##calculate the average KS d across samples
      ##maintain a variable with the average abundance in the selected window
      KS.matrix[idx,j] = mean(KS.vector[KS.vector != 0])
      abn.matrix[idx,j] = mean(sorted.matrix[idx:(idx+n.elements.per.window-1),j])
    }
  }
  returnObject <- list("abn" = abn.matrix, "KS.matrix" = KS.matrix)
  return(returnObject)
}

#' Function to calculate the Minkowski distance matrices using a sliding window
#'
#' This function generates an average Minkowski distance value for every sliding window
#' for each sample, in the expression matrix. 
#' @param expression.matrix expression matrix, normalized
#' @param p the Minkowski order; p=1 : Manhattan dist; p=2 : Euclidean dist; p=3 : Cubic dist
#' @keywords Minkowski distances on sliding windows
#' @export
#' @examples
calculateMinkowskiMatrix <- function(expression.matrix, p=2)
{
  cat("the input matrix has ",nrow(expression.matrix)," rows and ",ncol(expression.matrix)," cols\n")
  cat("number of genes: ",nrow(expression.matrix),"\n")
  cat("number of samples: ",ncol(expression.matrix),"\n")
  cat("the Minkowski distance of rank ",p," will be calculated\n")
  
  n.elements.per.window = floor(nrow(expression.matrix) / 10)
  cat("the number of elements per window is ",n.elements.per.window,"\n")
  
  Mnk.matrix = rep(0,ncol(expression.matrix)*(nrow(expression.matrix)-n.elements.per.window+1))
  Mnk.matrix = matrix(Mnk.matrix,ncol=ncol(expression.matrix))
  
  abn.matrix = rep(0,ncol(expression.matrix)*(nrow(expression.matrix)-n.elements.per.window+1))
  abn.matrix = matrix(abn.matrix,ncol=ncol(expression.matrix))
  
  for(j in 1:ncol(expression.matrix))
  {
    cat("working with sample ",j,"\n")
    sorted.matrix = expression.matrix[order(expression.matrix[,j]),]
    for(idx in 1:(nrow(expression.matrix)-n.elements.per.window+1))
    {
      #focus on a sliding window
      vector.Mnk = rep(0,ncol(expression.matrix))
      for(k in 1:ncol(expression.matrix))
      {
        if(j != k)
        {
          col.j = sorted.matrix[idx:(idx+n.elements.per.window-1),j]
          names(col.j)=rep("",length(col.j))
          col.k = sorted.matrix[idx:(idx+n.elements.per.window-1),k]
          names(col.k)=rep("",length(col.k))
          cor.in = cbind(col.j,col.k)
          if(p==1) {vector.Mnk[k] = dist(t(cor.in),method="manhattan")}
          if(p==2) {vector.Mnk[k] = dist(t(cor.in),method="euclidean")}
          if(p>2)  {vector.Mnk[k] = dist(t(cor.in),method="minkowski",p=p)}
        }#end if
      }#end for k
      Mnk.matrix[idx,j] = mean(vector.Mnk[vector.Mnk != 0])
      abn.matrix[idx,j] = mean(sorted.matrix[idx:(idx+n.elements.per.window-1),j])
    }#end for on the sliding windows
  }
  ##normalized on the abundance matrix
  Mnk.matrix = Mnk.matrix / abn.matrix
  returnObject <- list("abn" = abn.matrix, "dist.matrix" = Mnk.matrix)
  return(returnObject)
}

##plot all correlations
plot.PSK <- function(file.name, abn.matrix=NULL, P.matrix=NULL, S.matrix=NULL,K.matrix=NULL)
{
  if(is.null(abn.matrix)){cat("please include an abundance matrix\n");return;}
  if(is.null(P.matrix) && is.null(S.matrix) && is.null(K.matrix))
  {cat("please include at least one correlation matrix\n");return;}
  cat("number of samples to be plotted: ",ncol(abn.matrix),"\n")
  
  pdf(file.name,width=7,height=7)
  for(j in 1:ncol(abn.matrix))
  {
    if(!is.null(P.matrix))
    {
      plot(x=abn.matrix[,j],y=P.matrix[,j],type='l',ylim=c(-1,1),main=paste("Sample",j),
           xlab="Abn(linear)",ylab='Correlation Coefficient')
      par(new=TRUE)
    }
    
    if(!is.null(S.matrix))
    {
      plot(x=abn.matrix[,j],y=S.matrix[,j],type='l',ylim=c(-1,1),col='red',main=paste("Sample",j),
           xlab="Abn(linear)",ylab='Correlation Coefficient')
      par(new=TRUE)
    }
    
    if(!is.null(K.matrix))
    {
      plot(x=abn.matrix[,j],y=K.matrix[,j],type='l',ylim=c(-1,1),col='blue',main=paste("Sample",j),
           xlab="Abn(linear)",ylab='Correlation Coefficient')
    }
    for(l in seq(-0.5,1,0.25))
    {abline(h=l,lty=2)}
    
    legend("bottomright", legend=c("Pearson", "Spearman","Kendall"),
           fill=c("black","red", "blue"))
    par(new=FALSE)
  }
  dev.off()
}

##plot the smoothed correlation coefficient
##this can be applied to other distances too, think of a better name for the function
plot.CC.smooth <- function(file.name, abn.matrix = NULL, corr.matrix = NULL, min.y=-1,max.y=1)
{
  pdf(file.name,width=7,height=7)
  for(j in 1:ncol(abn.matrix))
  {
    loessMod10 <- loess(corr.matrix[,j] ~ abn.matrix[,j], span=0.10)
    smoothed10 <- predict(loessMod10)
    loessMod25 <- loess(corr.matrix[,j] ~ abn.matrix[,j], span=0.25)
    smoothed25 <- predict(loessMod25)
    loessMod50 <- loess(corr.matrix[,j] ~ abn.matrix[,j], span=0.50)
    smoothed50 <- predict(loessMod50)
    
    plot(x=abn.matrix[,j],y=corr.matrix[,j],type='l',ylim=c(min.y,max.y),main=paste("Sample",j),
         xlab="Abn(linear)",ylab='Distance')
    lines(smoothed10,x=abn.matrix[,j],col='red')
    lines(smoothed25,x=abn.matrix[,j],col='blue')
    lines(smoothed50,x=abn.matrix[,j],col='green')
    
    legend("bottomright", legend=c("Loess10", "Loess25","Loess50"),
           fill=c("red", "blue","green"))
  }
  dev.off()
}

##this function can be used for both the entropy tests and the Minkowski distances
plot.distance <- function(file.name, abn.matrix = NULL, dist.matrix = NULL)
{
  dist.matrix[is.na(dist.matrix)]=0
  max.overall = max(dist.matrix)
  pdf(file.name,width=7,height=7)
  for(j in 1:ncol(abn.matrix))
  {
    plot(x=abn.matrix[,j],y=dist.matrix[,j],type='l',ylim=c(0,max.overall),main=paste("Sample",j),
         xlab="Abn(linear)",ylab='Distance(linear)')
  }
  dev.off()
}

##variant of the previous function; allows the plotting of multiple distances
##intended for plotting the Mink distances
##or groupping the KL div, KS and a MInk
plot.distances <- function(file.name, abn.matrix = NULL, 
                           dist.matrix1 = NULL, d1.name = NULL, 
                           dist.matrix2 = NULL, d2.name = NULL,
                           dist.matrix3 = NULL, d3.name = NULL)
{
  max.overall = max(max(dist.matrix1), max(dist.matrix2), max(dist.matrix3))
  pdf(file.name,width=7,height=7)
  for(j in 1:ncol(abn.matrix))
  {
    if(!is.null(dist.matrix1))
    {
      plot(x=abn.matrix[,j],y=dist.matrix1[,j],type='l',ylim=c(0,max.overall),main=paste("Sample",j),
           xlab="Abn(linear)",ylab='Distance(linear)')
      par(new=TRUE)
    }
    if(!is.null(dist.matrix3))
    {
      plot(x=abn.matrix[,j],y=dist.matrix2[,j],type='l',ylim=c(0,max.overall),main=paste("Sample",j),
           xlab="Abn(linear)",ylab='Distance(linear)',col='red')
      par(new=TRUE)
    }
    if(!is.null(dist.matrix2))
    {
      plot(x=abn.matrix[,j],y=dist.matrix3[,j],type='l',ylim=c(0,max.overall),main=paste("Sample",j),
           xlab="Abn(linear)",ylab='Distance(linear)',col='blue')
      
    }
    legend("bottomright", legend=c(d1.name, d2.name,d3.name),
           fill=c("black","red", "blue"))
    par(new=FALSE)
  }
  dev.off()
}



##function to draw the PCA plots for a given matrix of expression
##the format is: each sample is presented on a row
##the abundances of genes are presented on columns 
##one column contains the abundances for a gene in all samples.
##first column is isolate name [,1]
##second column is phenotype   [,2]
plotPCA <- function(expression.matrix,title="")
{
  library("FactoMineR")
  library("factoextra")
  library("ggsci")
  
  expr.PCA <- prcomp(expression.matrix[,3:ncol(expression.matrix)],
                     center = TRUE, scale. = TRUE)
  
  staph.p <- fviz_pca_ind(expr.PCA,
                          geom.ind = "point", # show points only (nbut not "text")
                          col.ind = expression.matrix[,2], # color by groups
                          palette = c("#00AFBB", "#FC4E07"),
                          addEllipses = TRUE, # Concentration ellipses
                          ellipse.type = "confidence", # for confidence ellipses
                          legend.title = "Groups") 
  
  ggpubr::ggpar(staph.p,
                title = "Principal Component Analysis",
                subtitle = title,
                # caption = "Source: nanostring",
                xlab = "PC1", ylab = "PC2",
                legend.title = "Species", legend.position = "top",
                ggtheme = theme_gray())
}

##calculate the euclidean distance between centroids
calculateCentroidDist <- function(expression.matrix)
{
  expr.PCA <- prcomp(expression.matrix[,3:ncol(expression.matrix)],
                     center = TRUE, scale. = TRUE)
  PCs <- expr.PCA$x[,1:2] ##select the first two principal components
  
  centroids <- c(0,0)
  
  pheno <- as.factor(expression.matrix[,2])
  for(p in levels(pheno))
  {
    selected = PCs[expression.matrix[,2]==p,]
    pca.centroid <- apply(selected,2,mean)
    centroids = rbind(centroids,pca.centroid)
  }
  
  ##distance between centroids
  return.dist <- dist(centroids[2:3,],method="euclidean")
  
  pheno.dist = rep(-1,length(levels(pheno)))
  idx=1
  ##distance within a class
  for(p in levels(pheno))
  {
    selected = PCs[expression.matrix[,2]==p,]
    pca.centroid <- apply(selected,2,mean)
    avg.d <- rep(-1,nrow(selected))
    for(i in 1:nrow(selected))
    {
      avg.d[i] = dist(rbind(pca.centroid,selected[i,]),method="euclidean")
    }
    pheno.dist[idx] = mean(avg.d)
    idx = idx+1
  }
  
  returnObject <- list("centroid.dist" = return.dist, "withinGroup" = pheno.dist)
  return(returnObject)
}

##calculate the separation base don the PC1 and PC2
predictSampleOutlier_PCA <- function(expression.matrix)
{
  all.samples      <- calculateCentroidDist(expression.matrix)
  cat("centroid distance for all samples: ",all.samples$centroid.dist,
      "and within groups",all.samples$withinGroup,"\n")
  
  # restricted.input <- rep(-1,nrow(expression.matrix)) 
  for(i in 1:nrow(expression.matrix))
  {
    selected = expression.matrix[-i,]
    restricted.input = calculateCentroidDist(selected)
    
    cat(i,as.character(expression.matrix[i,1])," has centroid distance ",restricted.input$centroid.dist,
        "and within groups",restricted.input$withinGroup,"\n")
  }
}


##function to iterate through the samples and idnetify clear outliers
calculateCentroidDistances <- function(expression.matrix)
{
  pheno <- as.factor(expression.matrix[,2])
  for (p in levels(pheno))
  {
    selected <- expression.matrix[expression.matrix[,2]==p,]
    ##calculate centroid
    p.centroid <- apply(selected[,3:ncol(selected)],2,mean)
    ##calculate Minkowski/correlation distance
    for(i in 1:nrow(selected))
    {
      d.in <- rbind(selected[i,3:ncol(selected)],p.centroid)
      d.man <- dist(d.in,method="manhattan")
      d.euc <- dist(d.in,method="euclidean")
      d.cor <- cor(x=t(selected[i,3:ncol(selected)]),y=p.centroid,method="pearson")
      
      d.man = round(d.man,digits = 2)
      d.euc = round(d.euc,digits = 2)
      d.cor = round(d.cor,digits = 2)
      
      # cat(as.character(selected[i,1]),as.character(selected[i,2]),
      #     d.man,d.euc,d.cor,"\n",sep='\t')
       cat(as.character(selected[i,1]),as.character(selected[i,2]),
          d.cor,"\n",sep='\t')
    }
  }
}

##find the samples which might be outliers
##affect the centroid the most
predictSampleOutlier <- function(expression.matrix)
{
  pheno <- as.factor(expression.matrix[,2])
  
  ##the return data frame will consist of:
  ## original correlation with the centroid
  ## mean correlation for the remaining samples with the centroid (when the sample i is removed)
  
  cor.return.sample <- "" ##sample id
  cor.return.out    <- "" ##sample output
  cor.return.centr  <- rep(100,nrow(expression.matrix)) ## the correlations of the sample with the centroid
  cor.return.excl   <- rep(100,nrow(expression.matrix)) ## the mean correlation when the sample is excluded
  
  idx = 1
  for (p in levels(pheno))
  {
    selected <- expression.matrix[expression.matrix[,2]==p,]
    ##calculate centroid
    p.centroid <- apply(selected[,3:ncol(selected)],2,mean)
    ##calculate Minkowski/correlation distance
    for(i in 1:nrow(selected))
    {
      d.in <- rbind(selected[i,3:ncol(selected)],p.centroid)
      d.man <- dist(d.in,method="manhattan")
      d.euc <- dist(d.in,method="euclidean")
      d.cor <- cor(x=t(selected[i,3:ncol(selected)]),y=p.centroid,method="pearson")
      
      d.man = round(d.man,digits = 2)
      d.euc = round(d.euc,digits = 2)
      d.cor = round(d.cor,digits = 2)
      
      cor.return.sample[idx] = as.character(selected[i,1])
      cor.return.out[idx]    = as.character(selected[i,2])
      cor.return.centr[idx]  = d.cor
      
      exclusion.matrix <- selected[-i,]
      p.centroid.new   <-  apply(exclusion.matrix[,3:ncol(selected)],2,mean)
      all.corr         <- rep(100,nrow(exclusion.matrix))
      for(j in 1:nrow(exclusion.matrix))
      {
        d.cor <- cor(x=t(exclusion.matrix[j,3:ncol(exclusion.matrix)]),
                     y=p.centroid.new,method="pearson")
        all.corr[j] = d.cor
      }
      
      cor.return.excl[idx] = round(mean(all.corr),digits = 2)
      idx = idx + 1
    }
  }
  
  return.object <- data.frame(cor.return.sample,cor.return.out,cor.return.centr,cor.return.excl)
  return(return.object)
}