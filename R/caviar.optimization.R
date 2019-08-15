caviar.optimization <- function(r,date,q.level,empiricalQuantile,nInitialCond=10,nInitVec=1000,is.plot=FALSE){
  

  RQfval <- matrix(NA,nrow=nInitVec,ncol=1) 
  
  initialTargetVectors <- matrix(runif(nInitVec*3,min=0,max=1),nrow=nInitVec,ncol=3)
  
  for (i in 1:nInitVec){
    RQfval[i,] <- caviar.obj(initialTargetVectors[i, ], r, empiricalQuantile, q.level, 1)
  }
  
  Results          <- cbind(RQfval, initialTargetVectors)
  SortedResults    <- Results[order(Results[,1]),] # sort the results
  BestInitialCond  <- SortedResults[1:nInitialCond,-1] # keep the top nInitialCond
  est <- NULL
  for (dt in 1:nInitialCond){
    est[[dt]] <-  suppressWarnings(optimx::optimx(BestInitialCond[dt,],caviar.obj,r=r,empiricalQuantile=empiricalQuantile,q.level=q.level,out=1,method=c("Nelder-Mead")))
  }
  val <- NULL
  for(j in 1:length(est)){val[j] <- est[[j]]$value }
  temp <- est[which.min(val)] 
  coeff.CAViaR <- as.numeric(temp[[1]][1:dim(initialTargetVectors)[2]])
  rq.val.CAViaR <- as.numeric(temp[[1]][4])
  tmp <- caviar.obj(coeff.CAViaR,r,empiricalQuantile,q.level,2)
  cond.quant.CAViaR <- tmp$VaR
  if (is.plot){
    plot(date,r,type='l',main=paste0("CAViaR. Quantile level: (", q.level, ")"), xlab='Months',ylab='')
    lines(date,cond.quant.CAViaR,type='l',col="red")
  }
 return(list(coeff.CAViaR=coeff.CAViaR,rq.val.CAViaR=rq.val.CAViaR,cond.quant.CAViaR=cond.quant.CAViaR))
 }



caviar.obj <- function(beta,r,empiricalQuantile,q.level,out){
  
  # get VaR:
  VaR <- EvaluateVaR(beta,r,empiricalQuantile)
  
  # compute hit stat:
  Hit <- as.numeric(r < VaR) - q.level
  
  # compute regression quantile loss:
  RQ  <- - sum(Hit*(r - VaR))
  
  if ((RQ == Inf) || (RQ != RQ)){
    RQ <- 1e+100
  }
  
  if (out == 1){
    res <- RQ
  } else if (out==2){
    res <- list(VaR=VaR,Hit=Hit)
  }
  return(res)
}



EvaluateVaR <- function(beta,r,empiricalQuantile){

  lVaR <- length(r)
  VaR <- matrix(NA,nrow=lVaR,ncol=1)
  VaR[1,] <- empiricalQuantile

  # SAV CAViaR specification:

  for (i in 2:lVaR){
    VaR[i,] <- beta[1] + beta[2] * VaR[i - 1, ] + beta[3] * abs(r[i - 1])
  }

return(as.vector(VaR))
}