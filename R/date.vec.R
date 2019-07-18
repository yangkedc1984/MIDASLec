date.vec <-
function(s) {
  mat <- matrix(0, nrow = length(s), ncol = 6 )
  a <- as.POSIXlt(as.Date(s, "1970-01-01"))
  mat[,1] <- a$year + 1900
  mat[,2] <- a$mon + 1
  mat[,3] <- a$mday 
  mat[,4] <- a$hour
  mat[,5] <- a$min
  mat[,6] <- a$sec
  
  list(mat)
}
