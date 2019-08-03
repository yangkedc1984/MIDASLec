forecast.ta.example <- function(obj,data){
  
  out.y <- data$out.y
  out.x <- data$out.x
  data.oos <- cbind(rep(1,times=length(data$out.y)), data$out.lag.y, rowMeans(data$out.x[,1:3]),rowMeans(data$out.x[,4:6]),rowMeans(data$out.x[,7:9]))
  beta <- coefficients(obj)

  pred <- data.oos%*%beta
  rmse <- sqrt(mean((pred - out.y)^2))
  return(list(pred=pred,rmse=rmse))
}

