mode <-
function(data) {
  
  # mode: mode (most frequent value) of a vector
  #
  # Input Arguments:
  #
  # data: a vector of data
  
  # Output Arguments:
  
  # dataMode: mode of the vector
  #
  # countMax: frequency count at the mode
  
  
  nobs <- length(data)
  data <- sort(data)
  count <- 1
  countMax <- 1
  dataMode = data[1]
  
  for (t in 2:nobs){
    if(data[t]==data[t-1]){
      count <- count + 1
    }  
    if(data[t]!=data[t-1]){ 
      if(count > countMax){
        countMax <- count 
        dataMode <- data[t-1]
      } 
      count <- 1
    }
  } # end of for
  if (count > countMax) {
    countMax <- count
    dataMode  <- data[nobs]
  }
  return(list(dataMode = dataMode, countMax = countMax))
}
