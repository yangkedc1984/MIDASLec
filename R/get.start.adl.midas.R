get.start.adl.midas <- function (y, X, z = NULL, weight, par.num.weight, num.evals=1000, num.coef=10, is.intercept = TRUE) 
{
  d <- ncol(X)
  nw <- par.num.weight
  model <- na.omit(cbind(y, X, z))
  y <- model[, 1]
  XX <- model[, -1]
  n <- nrow(model)
  if (is.null(z)) {
    all_coef <- function(p) {
      weight(p, d)
    }
    p.eval <- matrix(NA,ncol=par.num.weight,nrow=num.evals)
    # get slope via OLS:
    p.eval[,1] <- rep(as.numeric(lm(y~rowMeans(X))$coef[2]),times=num.evals)+runif(num.evals,min=-0.5,max=0.5)
    p.eval[,2:par.num.weight] <- matrix(runif(num.evals*(par.num.weight-1),min=-0.5,max=0.5),ncol=par.num.weight-1,nrow=num.evals)
  } else {
    all_coef <- function(p) {
      c(weight(p[1:nw], d), p[-nw:-1])
    }
    p.eval <- matrix(NA,ncol=par.num.weight+1,nrow=num.evals)
    # get slope via OLS:
    p.eval[,1] <- rep(as.numeric(lm(y~rowMeans(X))$coef[2]),times=num.evals)+runif(num.evals,min=-0.5,max=0.5)
    p.eval[,2:par.num.weight] <- matrix(runif(num.evals*(par.num.weight-1),min=-0.5,max=0.5),ncol=par.num.weight-1,nrow=num.evals)
    p.eval[,par.num.weight+1] <- rep(as.numeric(lm(y~z)$coef[2]),times=num.evals)+runif(num.evals,min=-0.1,max=0.1)
    
    
    
  }
  fn0 <- function(p) {
    sum((y - XX %*% all_coef(p))^2)
  }
  f.eval <- NULL
  for (i in 1:num.evals){
    f.eval[i] <- fn0(p.eval[i,])
  }
  all <- cbind(f.eval,p.eval)
  all <- all[order(all[,1]),]
  coefs <- all[1:num.coef,-1]
  if (is.intercept){
    int <- rep(mean(y),times=num.coef)
    if (num.coef==1){
      coefs <- c(coefs,int)
    } else {
      coefs <- cbind(coefs,int)
    }
  }
  return(coefs)
}