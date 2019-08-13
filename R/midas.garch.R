midas.garch <- function(y,period,nlag,params0,regressor=NULL,rollingWindow=FALSE,thetaM=FALSE){
  
  # Replace missing values by the sample average
  y[is.na(y)] <- mean(y,na.rm = TRUE)
  Yfull <- y
  
  # Reshape the observation vector as a matrix, in which each column contains 
  # observations in a week/month/quarter/year
  nobs <- length(y)
  nMonth <- ceiling(nobs/period)
  Y <- matrix(NA,nrow=period,ncol=nMonth)
  
  Y[1:nobs] <- y 
  date.rv <- matrix(NA,nrow=period,ncol=nMonth)
  tmp <- nasdaq.trunc[,1]
  date.rv[1:nobs] <- tmp
  
  if (!is.null(regressor)){
    regressor[is.na(regressor)] <- mean(regressor,na.rm=TRUE)
    regressorFull <- regressor
    regressor <- regressor[1:nobs]
  }
  # Realized volatility is fixed in a week/month/quarter/year
  if (rollingWindow==TRUE){
    # Realized volatility varies every period due to rolling windows
    RV <- matrix(0,nrow=nobs,ncol=1) 
    for (t in period:nobs){
      RV[t] = sum(y[seq(t-period+1,t,by=1)]^2,na.rm = TRUE)
    }
    RV[1:(period-1)] <- RV[period]
  } else {
    if (is.null(regressor)){
      RV <- colSums(Y^2,na.rm = TRUE)
      # Odd days of the last week/month/quarter/year
      if (nobs %% period > 0) {
        RV[length(RV)] <- sum(Y[(length(Y)-period+1):length(Y)]^2,na.rm=TRUE)
      }
    } else {
      RegMat <- matrix(NA,nrow=period,ncol=nMonth) 
      RegMat[1:nobs] <- regressor
      RV <- colMeans(RegMat,na.rm=TRUE)
      # Odd days of the last week/month/quarter/year
      if (nobs %% period > 0) {
        RV[length(RV)] <- mean(regressor[seq(length(regressor)-period+1,length(regressor),by=1)],na.rm=TRUE) 
      }
    }
  }
  
  
  
  if (thetaM) {
    lb <- c(-1,0,0,-1,1.001,-1)
    ub <- c(1,1,1,1,50,1)
  } else  {
    lb <- c(-1,0,0,0,1.001,0,0)
    ub <- c(1,1,1,1,50,1)
  }
  
  myfun <- function(params,Y,RV,nlag,nobs,rollingWindow,thetaM){
    
    tmp <- fML(params,Y,RV,nlag,nobs,rollingWindow,thetaM)
    
    out <- -sum(tmp$logL)
    
    return(out)
  }
  
  est <- suppressWarnings(optimx::optimx(params0,myfun,Y=Y,RV=RV,nlag=nlag,nobs=nobs,rollingWindow=rollingWindow,thetaM=thetaM,method=c("L-BFGS-B"),lower=lb,upper=ub))
  if (est$value>1e100){
    est <- suppressWarnings(optimx::optimx(params0,myfun,Y=Y,RV=RV,nlag=nlag,nobs=nobs,rollingWindow=rollingWindow,thetaM=thetaM))
    idx <- which(est$value<1e100)
    if (length(idx)>1){idx <- 1 }
    est <- est[idx,]
  }
  estParams <- as.numeric(est[1:6])
  out <- fML(estParams,Y,RV,nlag,nobs,rollingWindow,thetaM)
  # AIC and BIC
  logLikeSum <- sum(out$logL,na.rm=TRUE)
  numParam <- length(params0)
  aic <- -2*logLikeSum + 2*numParam
  bic <- -2*logLikeSum + numParam*log(nobs)
  adjustSampleSize <- sum(out$logL!=0 * !is.na(out$logL))
  
  nobsBig <- length(Yfull)
  nMonthBig <- ceiling(nobsBig/period)
  Ybig <- matrix(NA,nrow=period,ncol=nMonthBig)
  Ybig[1:nobsBig] <- Yfull
  
  if (rollingWindow==TRUE){
    RVBig <- matrix(0,nrow=nobsBig,ncol=1)
    for (t in period:nobsBig){
      RVBig[t] <- sum(Yfull[seq(t-period+1,t,by=1)]^2,na.rm = TRUE) 
    }
    RVBig[1:(period-1)] <- RVBig[period]
  } else {
    if (is.null(regressor)){
      RVBig <- colSums(Ybig^2,na.rm = TRUE)
      # Odd days of the last week/month/quarter/year
      if (nobsBig %% period > 0) {
        RVBig[length(RVBig)] <- sum(Ybig[(length(Ybig)-period+1):length(Ybig)]^2,na.rm=TRUE)
      } 
    } else {
      RegMatBig <- matrix(NA,nrow=period,ncol=nMonth) 
      RegMatBig[1:nobsBig] <- regressorFull
      RVBig <- colMeans(RegMat,na.rm=TRUE)
      # Odd days of the last week/month/quarter/year
      if (nobsBig %% period > 0) {
        RVBig[length(RVBig)] <- mean(regressorFull[seq(length(regressorFull)-period+1,length(regressorFull),by=1)],na.rm=TRUE) 
      }
    }
  }
  out <- fML(estParams,Ybig,RVBig,nlag,nobs,rollingWindow,thetaM)
  # One-step-ahead in-sample forecast validation
  RealizedY2 <- (Yfull - estParams[1])^2
  forecastError <- out$Variance - RealizedY2
  estSampleRMSE = sqrt(mean(forecastError^2))
  
  return(list(estParams=estParams,logL=out$logL,Variance=out$Variance,LongRun=out$LongRun,estSampleRMSE=estSampleRMSE))
}

fML  <- function(params,Y,RV,nlag,nobs,rollingWindow,thetaM) {
  
  
  # Matrix dimension
  period <- dim(Y)[1]
  nMonth <- dim(Y)[2]
  
  
  # Allocate parameters
  mu0 <- params[1]
  alpha0 <- params[2]
  beta0 <- params[3]
  theta0 <- params[4]
  w10 <- params[5]
  w20 <- NULL
  m0 <- params[6]
  
  # GARCH positive constraint
  intercept <- 1 - alpha0 - beta0
  if ((intercept < 0) || (alpha0 < 0) || (beta0 < 0) || (alpha0 > 1) || (beta0 > 1)) {
    logL  <- matrix(Inf,nrow=nobs,ncol=1)
    Variance  <- matrix(NA,nrow=nobs,ncol=1)
    ShortRun  <- matrix(NA,nrow=nobs,ncol=1)
    LongRun  <- matrix(NA,nrow=nobs,ncol=1)
    return(list(logL=logL,Variance=Variance,ShortRun=ShortRun,LongRun=LongRun))
  }
  
  # theta and m are squared for compatibility with others' codes
  if (thetaM==FALSE){
    theta0 <- theta0 * theta0
    m0 <- m0 * m0
  }
  # Deflate observations and take squared residuals
  Ydeflate <- Y - mu0
  ResidSq <- Ydeflate * Ydeflate
  
  ShortRun <- matrix(1,nrow=period,ncol=nMonth)
  tauAvg <- m0 + theta0 * mean(RV,na.rm = TRUE)
  Variance <- tauAvg * matrix(1,nrow=period,ncol=nMonth)
  
  
  if (rollingWindow==TRUE){
    
    # Compute MIDAS weights
    nlagBig <- period*nlag
    weights <- midasBetaWeights(nlagBig,w10,w20)
    loopStart <- nlagBig + 1
    
    for (t in loopStart:nobs){
      
      # Compute long-run component
      # Refer to Eq (5) in Engle et al. (2013)
      tau <- m0 + theta0 * (sum(weights * RV[seq(t-nlagBig,t-1,by=1)],na.rm=TRUE))
      alphaTau <- alpha0 / tau
      
      # Compute short-run component
      # Refer to Eq (4) in Engle et al. (2013)
      ShortRun[t] <- intercept + alphaTau * ResidSq[t-1] + beta0 * ShortRun[t-1]
      # Compute conditional variance
      # Refer to Eq (3) in Engle et al. (2013)
      Variance[t] = tau * ShortRun[t]
      
    }
    
  } else if (rollingWindow==FALSE){
    
    weights <- t(midasBetaWeights(nlag,w10,w20))
    
    
    for (t in (nlag+1):nMonth){
      RVuse <- RV[seq(t-nlag,t-1,by=1)]
      
      # Compute long-run component
      # Refer to Eq (5) in Engle et al. (2013)
      tau <- m0 + theta0 * (sum(RVuse * weights))
      alphaTau <- alpha0 / tau
      
      # Compute short-run component
      # Refer to Eq (4) in Engle et al. (2013)
      for (n in 1:period) {
        ind <- (t-1)*period + n
        ShortRun[ind] <- intercept + alphaTau * ResidSq[ind-1] + beta0 * ShortRun[ind-1]
      }
      
      Variance[,t] <- tau * ShortRun[,t]
    }
  }
  
  if (any(Variance<0,na.rm=TRUE)){
    logL  <- matrix(Inf,nrow=nobs,ncol=1)
    Variance  <- matrix(NA,nrow=nobs,ncol=1)
    ShortRun  <- matrix(NA,nrow=nobs,ncol=1)
    LongRun  <- matrix(NA,nrow=nobs,ncol=1)
    return(list(logL=logL,Variance=Variance,ShortRun=ShortRun,LongRun=LongRun))
  }
  
  
  # Compute GARCH-MIDAS log likelihood 
  logLMatrix <- - 0.5 * ( log(2*pi*Variance) + ResidSq / Variance )
  logLMatrix[,1:nlag] <- 0
  
  logL <- logLMatrix[1:nobs]
  
  
  # If logMatrix are all NaNs, logL would be all zeros
  if ((all(logL == 0)) || sum(logL,na.rm=TRUE) == 0){
    logL = matrix(Inf,nrow=nobs,ncol=1)
  }
  
  Variance <- Variance[1:nobs]
  ShortRun <- ShortRun[1:nobs]
  LongRun <- Variance / ShortRun
  
  return(list(logL=logL,Variance=Variance,ShortRun=ShortRun,LongRun=LongRun))
}

midasBetaWeights <-function(nlag,param1,param2=NULL){
  seq. <- seq(nlag,1,by=-1)
  eps <- .Machine$double.eps
  if (is.null(param2)==TRUE){    
    weights <- (1-seq./nlag+10*eps)^(param1-1)
  } else {
    weights <- (1-seq./nlag+10*eps)^(param1-1) * (seq./nlag)^(param2-1)  
  }
  weights <- weights / sum(weights,na.rm=TRUE)
  return(weights)
}