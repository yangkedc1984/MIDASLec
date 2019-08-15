midas.optimization.rq <- function(dataY,dateY,optionsmidas,weights,q.level,nInitialCond=5,is.plot=FALSE) {
  # q.level - quantile level
  # nInitialCond - Number of initial conditions for the optimization set to 5

  
  N <- 1  # total # of regressors
  
  #______________________ Data construction _____________________________%
  #___________________________________________________________________%
  tmp <- datageneration.midas(dataY,optionsmidas,dateY)
  
  y <- tmp$y
  X <- abs(tmp$x)
  date <- tmp$dateselected
  #______________________ Initial guess _____________________________%
  #______________________________________________________________%
  # Select best starting points: select numInVec random initial values, calculate the obj fn at these points,
  # then select the best nInitialCond ones as candidate initial values from where to start optimizing
  
  set.seed(1)      # initialize the seed for replicability
  Ns <- 1000       # number of trial vectors
  
  
  if (weights=="beta"){
    # random vectors of MIDAS coeff
    initialTargetVectors <- matrix(NA,nrow=Ns,ncol=N*3)
    for (col in 1:N){ # for each regressor, if more than 1
      initialTargetVectors[,seq((col-1)*3+1,col*3,by=1)] = cbind(rnorm(Ns), matrix(runif(Ns*2,min=0,max=40),nrow=Ns,ncol=2))
    }
  } else if (weights=="betaconstr"){
    # random vectors of MIDAS coeff
    initialTargetVectors <- matrix(NA,nrow=Ns,ncol=N*2)
    for (col in 1:N){ # for each regressor, if more than 1
      initialTargetVectors[,seq((col-1)*2+1,col*2,by=1)] = cbind(rnorm(Ns), runif(Ns,min=0,max=40))
    }
    
  } else if (weights=="exp"){
    # random vectors of MIDAS coeff
    initialTargetVectors <- matrix(NA,nrow=Ns,ncol=N*3)
    for (col in 1:N){ # for each regressor, if more than 1
      initialTargetVectors[,seq((col-1)*3+1,col*3,by=1)] = cbind(rnorm(Ns), matrix(runif(Ns*2,min=-0.1,max=0),nrow=Ns,ncol=2))
    }
  }
  #  add the constant term
  initialTargetVectors <- cbind(rnorm(Ns)/10,initialTargetVectors)
  # value of obj fn at these points
  V <- matrix(NA,nrow=Ns,ncol=1)
  for (k in 1:Ns){
    V[k,1]  <-  midas.objfn.rq(initialTargetVectors[k,],y,X,1,weights,q.level)
  }          
  Results          <- cbind(V, initialTargetVectors)
  SortedResults    <- Results[order(Results[,1]),] # sort the results
  BestInitialCond  <- SortedResults[1:nInitialCond,-1] # keep the top nInitialCond
  
  #______________________ Actual estimation _____________________________%
  #___________________________________________________________________%
  
  est <- NULL
  for (dt in 1:nInitialCond){
    est[[dt]] <-  suppressWarnings(optimx::optimx(BestInitialCond[dt,],midas.objfn.rq,y=y,X=X,OUT=1,weights=weights,q.level=q.level,method=c("Nelder-Mead")))
  }
  val <- NULL
  for(j in 1:length(est)){val[j] <- est[[j]]$value }
  temp <- est[min(val)==val]
  coeff.MIDAS <- as.numeric(temp[[1]][1:dim(initialTargetVectors)[2]])
  rq.val.MIDAS <- min(val) 
  cond.quant.MIDAS <- midas.objfn.rq(coeff.MIDAS,y,X,2,weights,q.level)
  if (is.plot){
    plot(date,y,type='l',main=paste0("MIDAS quantile. Quantile level: (", q.level, ")"), xlab='Months',ylab='')
    lines(date,cond.quant.MIDAS,type='l',col="red")
  }
  
  return(list(coeff.MIDAS=coeff.MIDAS,rq.val.MIDAS=rq.val.MIDAS,cond.quant.MIDAS=cond.quant.MIDAS,y=y,date=date))
}


midas.objfn.rq <- function(theta,y,X,OUT,weights,q.level) {
  # objective fn for the MIDAS estimation
  # MIDAS vol
  val <- theta[1]+theta[2]*midas.filter(theta[-c(1,2)],X,weights)
  # rq function:
  hit.stat <- q.level-(y<val)
  rq.stat <- t(hit.stat)%*%(y-val)
  if ((rq.stat == Inf)  || (is.complex(rq.stat))) {
    rq.stat <- 1e+10000
  }
  if (OUT == 1){
    return(as.numeric(rq.stat))
  } else if (OUT == 2){
    return(val)
  }
}

midas.filter <- function(coeff,x,weights) {
  if (weights=="beta"){
    smpl <- dim(x)[2]
    be.values <- beta.weights(smpl,coeff[1],coeff[2])
    y.h <- x%*%be.values
    
  } else if  (weights=="betaconstr"){
    smpl <- dim(x)[2]
    be.values <- beta.weights(smpl,1,coeff[1])
    y.h <- x%*%be.values
    
  } else if (weights=="exp"){
    smpl <- dim(x)[2]
    be.values <- beta.weights(smpl,coeff[1],coeff[2])
    y.h <- x%*%be.values
    
  } 
  
  
  return(y.h)
}

beta.weights <- function(dayLag,k1,k2){
  eps <- .Machine$double.eps
  u <- seq(eps,1-eps,length.out = dayLag)
  
  k1 <- (k1*(k1>0)+1e-8*(!k1>0))*(k1<300)+300*(!k1<300) # 0<k1<300
  k2 <- (k2*(k2>0)+1e-8*(!k2>0))*(k2<300)+300*(!k2<300) # 0<k1<300
  
  be.values <- u^(k1-1)*(1-u)^(k2-1)
  be.values <- be.values/sum(be.values)
}


exp.weights <- function(dayLag,k1,k2){
  iii <- seq(1,dayLag,length.out = dayLag)
  poli <- exp(k1*iii+k2*iii^2)/sum(exp(k1*iii+k2*iii^2))
}


