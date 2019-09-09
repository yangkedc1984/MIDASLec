mixed.freq.data.single <-
function(data.refdate,data.x,data.xdate,x.lag,horizon,est.start,est.end,disp.flag=TRUE) {
  # complete data
  mask.na <- !is.na(data.refdate)
  data.refdate <- data.refdate[mask.na]
  mask.na <- !is.na(data.x)
  data.x <- data.x[mask.na]
  data.xdate <- data.xdate[mask.na]
  data.refdate <- as.Date(data.refdate)
  data.x <- as.vector(data.x)
  data.xdate <- as.Date(data.xdate)
  
  data.refdate.vec <- as.Date(data.refdate)
  data.xdate.vec <- as.Date(data.xdate)
  
  est.start <- as.Date(est.start)
  est.end <- as.Date(est.end)
  
  data.refdate.vec <- date.vec(data.refdate.vec)
  data.xdate.vec <- date.vec(data.xdate.vec)
  data.refdate.vec <- matrix(unlist(data.refdate.vec),nrow=length(data.refdate))
  data.xdate.vec <- matrix(unlist(data.xdate.vec),nrow=length(data.xdate))
  data.refdate.num <- data.refdate
  data.xdate.num <- data.xdate
  
  date.format = c('year(s)','month(s)','day(s)','hour(s)','minute(s)','second(s)')
  period.ref <- data.freq(data.refdate.vec)$period 
  unit.ref <- data.freq(data.refdate.vec)$unit 
  period.x <- data.freq(data.xdate.vec)$period  
  unit.x <- data.freq(data.xdate.vec)$unit
  
  ref.lag <- lag.num(1,period.ref,unit.ref)
  x.lag <- lag.num(x.lag,period.x,unit.x)
  horizon <- lag.num(horizon,period.x,unit.x)
  if (ref.lag < 0){
    stop('ref.lag cannot be negative.')
  }
  if (x.lag < 0) {
    stop('x.lag cannot be negative')
  }
  # Minimum and maximum dates that data support
  min.date.ref <- data.refdate.num[ref.lag+1]
  min.date.x <- data.xdate.num[max(1,x.lag+horizon)]
  if (min.date.ref > min.date.x){
    min.date <- min.date.ref
  } else {
    min.date <- min.date.x
  }
  max.date.ref <- data.refdate.num[length(data.refdate.num)]
  max.date.x = data.xdate.num[length(data.xdate.num)]
  if (horizon < 0){
    max.date.x <- data.xdate.vec[dim(data.xdate.vec)[1],]
    max.date.x[unit.x] <- max.date.x[unit.x] + period.x * horizon
    max.date.x <- ISOdate(max.date.x[1],max.date.x[2],max.date.x[3],max.date.x[4],max.date.x[5],max.date.x[6])
    max.date.x <- as.Date(max.date.x)
  }
  if (max.date.ref > max.date.x){
    max.date <- max.date.x
  } else {
    max.date <- max.date.ref
  }
  # Check and set default sample period
  if (is.null(est.start)){
    est.start <- min.date
  } else { if(est.start < min.date) {#warning('Start date cannot be earlier than possible due to lagged regressors. Reset start date to most recent possible: ', paste(min.date))
    est.start <- min.date}
  }
  if (is.null(est.end)){
    est.end <- max.date
  } else { if(est.end > max.date) {warning('Terminal date cannot be later than largest date accounting for lags. Reset to largest date possible: ', paste(max.date))
    est.end <- max.date}
  }
  # Construct reference date data
  tol <- 1e-10
  loc.start <- min(which((data.refdate.num >= est.start-tol) == TRUE))
  loc.end <- min(which((data.refdate.num >= est.end-tol) == TRUE)) 
  est.refdate <- data.refdate.num[loc.start:loc.end]
  
  loc.forecast.end <- min(which((data.refdate.num >= max.date-tol) == TRUE))
  if(loc.end+1<=loc.forecast.end){
    out.refdate <- data.refdate.num[seq(loc.end+1,loc.forecast.end,by=1)]
    n.forecast <- length(out.refdate)
  } else {
    out.refdate <- NULL
    n.forecast <- length(out.refdate)
  }
  nobs <- loc.end - loc.start + 1
  
  est.x <- est.xdate <- matrix(NaN,nrow=nobs,ncol=x.lag) 
  for (t in 1:nobs){
    loc <- min(which((data.xdate.num >= est.refdate[t]-tol) == TRUE)) 
    if (is.null(loc)) {
      loc <- length(data.xdate.num)
    }
    
    if(loc-horizon > length(data.x)){    
      nobs <- t - 1
      #est.ref = est.ref[seq(1,nobs,1)]
      est.refdate = est.refdate[seq(1,nobs,1)]
      #est.lag.ref = est.lag.ref[seq(1,nobs,1)]
      #est.lag.refdate = est.lag.refdate[seq(1,nobs,1)]
      est.x = est.x[seq(1,nobs,1)]
      est.xdate = est.xdate[seq(1,nobs,1)]
      max.date = est.refdate[length(est.refdate)]
      warning('Horizon is a large negative number. Observations are further truncated to max date possible: ', paste(max.date))
      break
    } else  {      
      est.x[t,] <- data.x[seq(loc-horizon,loc-horizon-x.lag+1,-1)]
      est.xdate[t,] <- data.xdate.num[seq(loc-horizon,loc-horizon-x.lag+1,-1)]
    }
  }
  if(loc.end+1<=loc.forecast.end){
    out.x <- out.xdate <- matrix(NaN,nrow=n.forecast,ncol=x.lag) 
    for(t in 1:n.forecast){
      loc <- min(which((data.xdate.num >= out.refdate[t]-tol) == TRUE))  
      if (is.null(loc)) {
        loc <- length(data.xdate.num)
      }
      
      if(loc-horizon > length(data.x)){      
        n.forecast <- t - 1
        #out.ref = out.ref[seq(1,n.forecast,1)]
        out.refdate = out.refdate[seq(1,n.forecast,1)]
        #out.lag.ref = out.lag.ref[seq(1,n.forecast,1)]
        out.lag.refdate = out.lag.refdate[seq(1,n.forecast,1)]
        out.x = out.x[seq(1,n.forecast,1)]
        out.xdate = out.xdate[seq(1,n.forecast,1)]
        break
      } else {
        out.x[t,] <- data.x[seq(loc-horizon,loc-horizon-x.lag+1,-1)] 
        out.xdate[t,] <- data.xdate.num[seq(loc-horizon,loc-horizon-x.lag+1,-1)] 
      } 
    }
  } else {
    out.x <- out.xdate <- NULL
  }
  
  if (disp.flag==T){
    # Display mixed frequency data
    cat('Frequency of Reference Date:',period.ref,date.format[unit.ref], "\n") 
    cat('Frequency of Data X:',period.x,date.format[unit.x], "\n") 
    cat('Start Date: ', paste(est.start), "\n") 
    cat('Terminal Date: ', paste(est.end), "\n") 
    
    # Display timeframe of mixed frequency regression
    cat('Mixed frequency data structure time frame:', "\n") 
    for(m in c(1,2,nobs)){
      cat("\n")
      cat(paste('Ref date(',as.Date(est.refdate[m],origin="1970-01-01"),')`s on: ', sep="")) 
      if (x.lag == 1) {
        cat(paste(' X(',as.Date(est.xdate[m],origin="1970-01-01"),')`s', sep="")) 
      }
      if (x.lag == 2){
        cat(paste(' X(',as.Date(est.xdate[m,1],origin="1970-01-01"),')`s X(',as.Date(est.xdate[m,dim(est.xdate)[2]],origin="1970-01-01"),')`s', sep=""))  
      }
      if (x.lag == 3) {
        cat(paste(' X(',as.Date(est.xdate[m,1],origin="1970-01-01"),')`s X(',as.Date(est.xdate[m,2],origin="1970-01-01"),')`s X(',as.Date(est.xdate[m,dim(est.xdate)[2]],origin="1970-01-01"),')`s', sep=""))  
      }
      if (x.lag > 3){
        cat(paste(' X(',as.Date(est.xdate[m,1],origin="1970-01-01"),')`s X(',as.Date(est.xdate[m,2],origin="1970-01-01"),')`s ... X(',as.Date(est.xdate[m,dim(est.xdate)[2]],origin="1970-01-01"),')`s', sep=""))  
      }
    }
  }
  output = list(est.refdate = est.refdate, est.x = est.x, est.xdate = est.xdate,
                out.refdate = out.refdate, out.x = out.x, out.xdate = out.xdate,
                x.lag = x.lag,min.date = min.date, max.date = max.date)
  return(output)
}
