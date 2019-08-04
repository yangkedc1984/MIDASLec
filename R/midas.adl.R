midas.adl <- function(data.y, data.ydate, data.x, data.xdate, est.start,est.end, horizon, x.lag, y.lag, polynomial, method, disp.flag,num.evals=10000,num.coef=5){
  #   polynomial  functional form of weights. Its value could be
  #                 o "nbeta": Beta polynomial
  #                 o "nealmon": Exp Almon polynomial
  #                 o "umidas":   Polynomial with step functions (Umidas)
  #                 o "timeaverage": time-averaged high-frequency data 
  #    methods for forecasting, method: 
  #                 o "fixed": fixed scheme 
  #                 o "rolling": rolling window
  #                 o "expand": expanding window
  # all other input parameters are as in MATALB.
  mf.data <- mixed.freq.data(data.y,data.ydate,data.x,data.xdate,x.lag,y.lag,horizon,est.start,est.end,disp.flag)
  
  est.y <- mf.data$est.y
  est.x <- mf.data$est.x
  est.lag.y <- mf.data$est.lag.y
  
  est.ydate <- mf.data$est.ydate
  est.xdate <- mf.data$est.xdate
  
  out.y <- mf.data$out.y
  out.x <- mf.data$out.x
  out.ydate <- mf.data$out.ydate
  out.xdate <- mf.data$out.xdate
  out.lag.y <- mf.data$out.lag.y


  nobs <- length(est.y)
  nforecast <- length(out.y)
  
  # methods:
  
  if (method=="fixed"){
    est.obj <- midas.estimate(est.y,est.x,est.lag.y,est.xdate,polynomial,num.evals,num.coef)
    pred.obj <- midas.forecast(est.obj,out.y,out.x,out.lag.y,polynomial)
  } else {
    nroll <- nforecast
    if (nroll == 0){
      stop('Rolling window does not apply because there are no rolling periods. Decrease "EstEnd".')
    }
    y.big <- c(est.y,out.y)
    x.big <- rbind(est.x,out.x)
    lag.y.big <-  rbind(est.lag.y,out.lag.y)
    x.date.big <- rbind(est.xdate,out.xdate)
    y.date.big <- c(est.ydate,out.ydate)
    pred <- matrix(0,nrow=nroll,ncol=1)
    for (t in 1:nroll){
      if (method=="rolling"){
        est.y.roll <- y.big[t:nobs-1+t]
        est.x.roll <- x.big[t:nobs-1+t,]
        est.lag.y.roll <- lag.y.big[t:nobs-1+t,]
        est.date.roll <- x.date.big[t:nobs-1+t,]             
      } else { 
        if (method=="expand"){
          est.y.roll <- y.big[1:nobs-1+t]
        est.x.roll <- x.big[1:nobs-1+t,]
        est.lag.y.roll <- lag.y.big[1:nobs-1+t,]
        est.x.date.roll <- x.date.big[1:nobs-1+t,]    
        } else {
          stop('method should be set to either: fixed, rolling, expand. Check!')
        }
      }
      out.y.roll <- y.big[nobs+t]
      out.x.roll <- x.big[nobs+t,]
      out.lag.y.roll <- lag.y.big[nobs+t,]
      out.y.dateroll <- y.date.big[nobs+t]
      if (t == 1){
        est.obj <- midas.estimate(est.y.roll,est.x.roll,est.lag.y.roll,est.x.date.roll,polynomial,num.evals,num.coef)
      } else  {
        est.obj <- midas.estimate(est.y.roll,est.x.roll,est.lag.y.roll,est.x.date.roll,polynomial,num.evals,num.coef=1,est.obj$coefficients)
      }
      tmp <- midas.forecast(est.obj,out.y.roll,t(as.matrix(out.x.roll)),t(as.matrix(out.lag.y.roll)),polynomial)
      pred[t] <- tmp$pred
    }
    pred.obj <- NULL
    pred.obj$pred <- pred
    pred.obj$rmse <- sqrt(mean((out.y-pred)^2))
  }
  return(list(est.obj=est.obj,pred.obj=pred.obj))
}

midas.estimate <- function(est.y,est.x,est.lag.y,est.xdate,polynomial,num.evals=10000,num.coef=5,startx.all=NULL){
  # nls-midas:
  if (polynomial=="nbeta") {
    weight <- nbeta
    if(is.null(startx.all)){# generate initial param guess:
    set.seed(123)
    startx.all <- get.start.adl.midas(y=est.y,X=est.x,z=est.lag.y,weight=weight,par.num.weight=3,num.evals=num.evals,num.coef=num.coef)
    }
    if (num.coef>1){
    est <- vals <- NULL
     for (j in 1:dim(startx.all)[1]) { 
       est[[j]] <- midas_r_plain(y=est.y,X=est.x,z=cbind(est.lag.y,rep(1,times=length(est.y))),weight=weight,startx=startx.all[j,-((dim(startx.all)[2]-1):dim(startx.all)[2])],startz=startx.all[j,((dim(startx.all)[2]-1):dim(startx.all)[2])],control=list(maxit=500))
       vals[j] <- est[[j]]$opt$value[1]
      }
      est.obj <- est[[which(min(vals)==vals)]]
    } else {
      est.obj <- midas_r_plain(y=est.y,X=est.x,z=cbind(est.lag.y,rep(1,times=length(est.y))),weight=weight,startx=startx.all[-c(4,5)],startz=startx.all[c(4,5)],control=list(maxit=500))
    }
    
  } else if (polynomial=="nealmon") {
    weight <- nealmon
    if(is.null(startx.all)){# generate initial param guess:
      set.seed(123)
      startx.all <- get.start.adl.midas(y=est.y,X=est.x,z=est.lag.y,weight=weight,par.num.weight=3,num.evals=num.evals,num.coef=num.coef)
    } 
    if (num.coef>1){
      est <- vals <- NULL
      for (j in 1:dim(startx.all)[1]) { 
        est[[j]] <- midas_r_plain(y=est.y,X=est.x,z=cbind(est.lag.y,rep(1,times=length(est.y))),weight=weight,startx=startx.all[j,-((dim(startx.all)[2]-1):dim(startx.all)[2])],startz=startx.all[j,((dim(startx.all)[2]-1):dim(startx.all)[2])],control=list(maxit=500))
        vals[j] <- est[[j]]$opt$value[1]
      }
      est.obj <- est[[which(min(vals)==vals)]]
    } else {
      est.obj <- midas_r_plain(y=est.y,X=est.x,z=cbind(est.lag.y,rep(1,times=length(est.y))),weight=weight,startx=startx.all[-c(4,5)],startz=startx.all[c(4,5)],control=list(maxit=500))
    }
  } else if (polynomial=="timeaverage"){
    est.obj <- lm(est.y~est.lag.y+rowMeans(est.x[,1:3])+rowMeans(est.x[,4:6])+rowMeans(est.x[,7:9]))
  } else if (polynomial=="umidas") {
    est.obj <- lm(est.y~est.lag.y+est.x)
  }
  return(est.obj)
}

midas.forecast <- function(obj,out.y,out.x,out.lag.y,polynomial){
  data <- NULL
  data$out.y <- out.y
  data$out.x <- out.x
  data$out.lag.y <- out.lag.y
  if (polynomial=="nbeta") {
    weight <- nbeta
    out <- forecast.adl(obj,weight=weight,par.num.weight=3,data,is.intercept=TRUE)
  } else if (polynomial=="nealmon") {
    weight <- nealmon
    out <- forecast.adl(obj,weight=weight,par.num.weight=3,data,is.intercept=TRUE)
  } else if (polynomial=="timeaverage"){
    out <- forecast.ta.example(obj,data)
  } else if (polynomial=="umidas") {
    out <- forecast.umidas(obj,data,is.intercept=TRUE)
  }
  return(out)
}
