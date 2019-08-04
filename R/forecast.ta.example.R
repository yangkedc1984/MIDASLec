forecast.ta.example <- function(obj,data){
  
  out.y <- data$out.y
  out.x <- data$out.x
  if (is.vector(out.x)){
    out.x <- t(as.matrix(out.x))
  }
  x <- cbind(rep(1,times=length(data$out.y)), data$out.lag.y, rowMeans(data$out.x[,1:3]),rowMeans(data$out.x[,4:6]),rowMeans(data$out.x[,7:9]))
  beta <- coefficients(obj)

  pred <- x%*%beta
  rmse <- sqrt(mean((pred - out.y)^2))
  return(list(pred=pred,rmse=rmse))
}

