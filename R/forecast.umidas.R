forecast.umidas <- function(obj,data,is.intercept=TRUE) {
  out.y <- data$out.y
  out.x <- data$out.x
  beta <- coefficients(obj)
  if(!is.null(data$out.lag.y)){
    x <- data$out.lag.y
  } else {
    x <- NULL
  }
  x <- cbind(x,data$out.x)
  if (is.intercept){
    x <- cbind(rep(1,times=length(out.y)) , x)
  }
  pred <- x%*%beta
  rmse <- sqrt(mean((pred - out.y)^2))
  return(list(pred=pred,rmse=rmse))
}