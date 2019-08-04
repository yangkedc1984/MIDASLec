forecast.adl <- function(obj,weight,par.num.weight,data,is.intercept=TRUE){
  
  out.y <- data$out.y
  out.x <- data$out.x
  
  beta <- coefficients(obj)
  # weight coefficients
  beta.w <- beta[1:par.num.weight]
  W <- weight(beta.w, d=dim(out.x)[2])
  
  if (is.intercept){
    c <- beta[length(beta)]
    beta <- beta[-length(beta)]
  }
  if (!is.null(data$out.lag.y)){
    yhl <- beta[(length(beta)-dim(data$out.lag.y)[2]+1):length(beta)]*data$out.lag.y
  } else {
    yhl <- 0
  }
  xhl <- out.x%*%W
  pred <- rep(c,times=length(out.y)) + xhl + yhl
  rmse <- sqrt(mean((pred - out.y)^2))
  return(list(pred=pred,rmse=rmse))
}
  
  