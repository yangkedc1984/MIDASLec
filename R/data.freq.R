data.freq <-
function(DateVec) {
  # data.freq: Identify data frequency
  #
  # Input Arguments:
  #  %
  # DateVec: T-by-6 R vector format data: [year,month,day,hour,min,sec]
  #
  # Output Arguments:
  #  %
  #% period: length of two consecutive dates
  #%
  #% unit: unit of length measure
  #%       o 1 = year
  #%       o 2 = month
  #%       o 3 = day
  #%       o 4 = hour
  #%       o 5 = minutes
  #%       o 6 = seconds
  #%
  #% Notes:
  #  %
  #% Frequency   period   unit
  #% yearly         1      1  
  #% semiannual     6      2 
  #% quarterly      3      2 
  #% monthly        1      2 
  #% biweekly       14     3
  #% weekly         7      3
  #% daily          1      3
  #% hourly         1      4
  #% minutely       1      5
  #% secondly       1      6
  
  DateDiff <- as.matrix(diff(DateVec))
  
  # Check annual or lower frequency
  modeUse = mode(DateDiff[,1])$dataMode
  if(modeUse >= 1) {
    period <- modeUse
    unit <- 1
    return(list(period = period,unit = unit))
  }
  
  # Check monthly frequency, quarter = 3 months, semiannual = 6 months
  modeUse <- mode(DateDiff[,2])$dataMode
  mask <- isTRUE(modeUse < 0)
  modeUse[mask] <- modeUse[mask] + 12
  if(modeUse >= 1){
    period <- modeUse
    unit <- 2
    return(list(period = period,unit = unit))
  }
  
  # Check daily frequency, week = 7 days, biweekly = 14 days
  modeUse <- mode(DateDiff[,3])$dataMode
  mask <- isTRUE(modeUse < 0)
  modeUse[mask] <- modeUse[mask] +30
  if(modeUse >= 1){
    period <- modeUse
    unit <- 3
    return(list(period = period,unit = unit))
  }
  
  # Check hourly frequency
  
  modeUse <- mode(DateDiff[,4])$dataMode
  mask <- isTRUE(modeUse < 0)
  modeUse[mask] <- modeUse[mask] + 24
  if(modeUse >= 1){
    period <- modeUse
    unit <- 4
    return(list(period = period,unit = unit))
  }
  
  # Check minutely frequency
  modeUse <- mode(DateDiff[,5])$dataMode
  mask <- isTRUE(modeUse < 0)
  modeUse[mask] <- modeUse[mask] + 60
  if(modeUse >= 1){
    period <- modeUse
    unit <- 5
    return(list(period = period,unit = unit))
  }
  
  
  # Check secondly frequency
  elapse <- diff.time.mf(DateVec[2:dim(DateVec)[1],6],DateVec[1:dim(DateVec)[1]-1,6],origin = "1970-01-01",units = "secs")
  period <- mean(elapse)
  unit   <- 6
  
  
  
}
