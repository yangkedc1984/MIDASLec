lag.num <-
function(x.lag,period,unit){
  
  if(is.numeric(x.lag)==T && is.atomic(x.lag)==T && length(x.lag) == 1L) {
    return(x.lag)
  }
  multiplier <- as.double(substr(x.lag, start = 1, stop = nchar(x.lag)-1))
  if (is.na(multiplier)) {
    stop('The description of lags cannot be recognized. The format should be 3m, 1q, etc')
  }
  #Convert multiplier to daily frequency (business days)
  ndaysPerYear <- 264
  ndaysPerQuarter <- 66
  ndaysPerMonth <- 22
  nhoursPerDay <- 8
  if(substr(x.lag, start = nchar(x.lag), stop = nchar(x.lag)) == 'y'){
    multiplier <- multiplier * ndaysPerYear
  }
  if(substr(x.lag, start = nchar(x.lag), stop = nchar(x.lag)) == 'q'){
    multiplier <- multiplier * ndaysPerQuarter
  }
  if(substr(x.lag, start = nchar(x.lag), stop = nchar(x.lag)) == 'm'){
    multiplier <- multiplier * ndaysPerMonth
  }
  if(substr(x.lag, start = nchar(x.lag), stop = nchar(x.lag)) == 'd'){
    multiplier <- multiplier * 1
  }
  if(substr(x.lag, start = nchar(x.lag), stop = nchar(x.lag)) == 'h'){
    multiplier <- multiplier / nhoursPerDay
  }
  if(substr(x.lag, start = nchar(x.lag), stop = nchar(x.lag)) == 's'){
    multiplier <- multiplier /  (nhoursPerDay*60*60)
  }
  if(unit == 1) {
    x.lag <- round(multiplier / (ndaysPerYear * period))
  }
  if(unit == 2) {
    x.lag <- round(multiplier / (ndaysPerMonth * period))
  }
  if(unit == 3) {
    x.lag <- round(multiplier / period)
  }
  if(unit == 4) {
    x.lag <- round(multiplier / (period / nhoursPerDay))
  }
  if(unit == 5) {
    x.lag <- round(multiplier / (period / nhoursPerDay / 60))
  }
  if(unit == 6) {
    x.lag <- round(multiplier / (period / nhoursPerDay / 60 / 60))
  }
  return(x.lag)
}
