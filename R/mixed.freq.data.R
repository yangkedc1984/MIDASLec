mixed.freq.data <-
function(data.y,data.ydate,data.x,data.xdate,x.lag,y.lag,horizon,est.start,est.end,disp.flag=TRUE) {
  # complete data
  mask.na <- !is.na(data.y)
  data.y <- data.y[mask.na]
  data.ydate <- data.ydate[mask.na]
  mask.na <- !is.na(data.x)
  data.x <- data.x[mask.na]
  data.xdate <- data.xdate[mask.na]
  data.y <- as.vector(data.y)
  data.ydate <- as.Date(data.ydate)
  data.x <- as.vector(data.x)
  data.xdate <- as.Date(data.xdate)
  
  
  data.ydate.vec <- as.Date(data.ydate)
  data.xdate.vec <- as.Date(data.xdate)
  
  
  est.start <- as.Date(est.start)
  est.end <- as.Date(est.end)
  
  
  data.ydate.vec <- date.vec(data.ydate.vec)
  data.xdate.vec <- date.vec(data.xdate.vec)
  data.ydate.vec <- matrix(unlist(data.ydate.vec),nrow=length(data.ydate))
  data.xdate.vec <- matrix(unlist(data.xdate.vec),nrow=length(data.xdate))
  data.ydate.num <- data.ydate
  data.xdate.num <- data.xdate
  
  date.format = c('year(s)','month(s)','day(s)','hour(s)','minute(s)','second(s)')
  period.y <- data.freq(data.ydate.vec)$period 
  unit.y <- data.freq(data.ydate.vec)$unit 
  period.x <- data.freq(data.xdate.vec)$period  
  unit.x <- data.freq(data.xdate.vec)$unit
  
  
  y.lag <- lag.num(y.lag,period.y,unit.y)
  x.lag <- lag.num(x.lag,period.x,unit.x)
  horizon <- lag.num(horizon,period.x,unit.x)
  if (y.lag < 0){
    stop('y.lag cannot be negative.')
  }
  if (x.lag < 0) {
    stop('x.lag cannot be negative')
  }
  
  # Minimum and maximum dates that data support
  min.date.y <- data.ydate.num[y.lag+1]
  min.date.x <- data.xdate.num[max(1,x.lag+horizon)]
  if (min.date.y > min.date.x){
    min.date <- min.date.y
  } else {
    min.date <- min.date.x
  }
  max.date.y <- data.ydate.num[length(data.ydate.num)]
  max.date.x = data.xdate.num[length(data.xdate.num)]
  if (horizon < 0){
    max.date.x <- data.xdate.vec[dim(data.xdate.vec)[1],]
    max.date.x[unit.x] <- max.date.x[unit.x] + period.x * horizon
    max.date.x <- ISOdate(max.date.x[1],max.date.x[2],max.date.x[3],max.date.x[4],max.date.x[5],max.date.x[6])
    max.date.x <- as.Date(max.date.x)
  }
  if (max.date.y > max.date.x){
    max.date <- max.date.x
  } else {
    max.date <- max.date.y
  }
  # Check and set default sample period
  if (is.null(est.start)){
    est.start <- min.date
  } else { if(est.start < min.date) {warning('Start date cannot be earlier than possible due to lagged regressors. Reset start date to most recent possible.')
    est.start <- min.date}
  }
  if (is.null(est.end)){
    est.end <- max.date
  } else { if(est.end > max.date) {warning('Terminal date cannot be later than largest date account for lags. Reset to largest date.')
    est.end <- max.date}
  }
  # Construct Y data
  tol <- 1e-10
  loc.start <- min(which((data.ydate.num >= est.start-tol) == TRUE))
  loc.end <- min(which((data.ydate.num >= est.end-tol) == TRUE)) 
  est.y <- data.y[loc.start:loc.end]
  est.ydate <- data.ydate.num[loc.start:loc.end]
  
  loc.forecast.end <- min(which((data.ydate.num >= max.date-tol) == TRUE))
  if(loc.end+1<=loc.forecast.end){
    out.y <- data.y[seq(loc.end+1,loc.forecast.end,by=1)]
    out.ydate <- data.ydate.num[seq(loc.end+1,loc.forecast.end,by=1)]
    n.forecast <- length(out.y)
  } else {
    out.y <- out.ydate <- NULL
    n.forecast <- length(out.y)
  }
  nobs <- loc.end - loc.start + 1
  # Construct lagged Y data
  est.lag.y <- est.lag.ydate <- matrix(NaN,nrow=nobs,ncol=y.lag)
  for (m in 1:y.lag){
    est.lag.y[,m] <- data.y[seq(loc.start-m,loc.end-m,1)]
    est.lag.ydate[,m] <- data.ydate.num[seq(loc.start-m,loc.end-m,1)]
  }
  
  if(loc.end+1<=loc.forecast.end){
    out.lag.y <- out.lag.ydate <- matrix(NaN,nrow=n.forecast,ncol=y.lag) 
    for (m in 1:y.lag){
      out.lag.y[,m] <- data.y[seq(loc.end-m,loc.forecast.end-m,1)]
      out.lag.ydate[,m] <- data.ydate.num[seq(loc.end-m,loc.forecast.end-m,1)]
    }
  } else {
    out.lag.y <- out.lag.ydate <- NULL
  }
  
  est.x <- est.xdate <- matrix(NaN,nrow=nobs,ncol=x.lag) 
  for (t in 1:nobs){
    loc <- min(which((data.xdate.num >= est.ydate[t]-tol) == TRUE)) 
    if (is.null(loc)) {
      loc <- length(data.xdate.num)
    }
    
    if(loc-horizon > length(data.x)){    
      nobs <- t - 1
      est.y = est.y[seq(1,nobs,1)]
      est.ydate = est.ydate[seq(1,nobs,1)]
      est.lag.y = est.lag.y[seq(1,nobs,1)]
      est.lag.ydate = est.lag.ydate[seq(1,nobs,1)]
      est.x = est.x[seq(1,nobs,1)]
      est.xdate = est.xdate[seq(1,nobs,1)]
      max.date = est.ydate[length(est.ydate)]
      warning('Horizon is a large negative number. Observations are further truncated to max date possible')
      break
    } else  {      
      est.x[t,] <- data.x[seq(loc-horizon,loc-horizon-x.lag+1,-1)]
      est.xdate[t,] <- data.xdate.num[seq(loc-horizon,loc-horizon-x.lag+1,-1)]
    }
  }
  if(loc.end+1<=loc.forecast.end){
    out.x <- out.xdate <- matrix(NaN,nrow=n.forecast,ncol=x.lag) 
    for(t in 1:n.forecast){
      loc <- min(which((data.xdate.num >= out.ydate[t]-tol) == TRUE))  
      if (is.null(loc)) {
        loc <- length(data.xdate.num)
      }
      
      if(loc-horizon > length(data.x)){      
        n.forecast <- t - 1
        out.y = out.y[seq(1,n.forecast,1)]
        out.ydate = out.ydate[seq(1,n.forecast,1)]
        out.lag.y = out.lag.y[seq(1,n.forecast,1)]
        out.lag.ydate = out.lag.ydate[seq(1,n.forecast,1)]
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
    cat('Frequency of Data Y:',period.y,date.format[unit.y], "\n") 
    cat('Frequency of Data X:',period.x,date.format[unit.x], "\n") 
    cat('Start Date: ', paste(est.start), "\n") 
    cat('Terminal Date: ', paste(est.end), "\n") 
    
    # Display timeframe of mixed frequency regression
    cat('Mixed frequency regression time frame:', "\n") 
    for(m in c(1,2,nobs)){
      cat("\n")
      cat(paste('Reg Y(',as.Date(est.ydate[m],origin="1970-01-01"),')`s on: ', sep="")) 
      if (y.lag == 1){
        cat(paste('Y(',as.Date(est.lag.ydate[m,],origin="1970-01-01"),')`s', sep="")) 
      }
      if (y.lag == 2){
        cat(paste('Y(',as.Date(est.lag.ydate[m,],origin="1970-01-01"),')`s Y(',as.Date(est.lag.ydate[m,dim(est.lag.ydate)[2]],origin="1970-01-01"),')`s', sep=""))  
      }
      if (y.lag > 2){
        cat(paste('Y(',as.Date(est.lag.ydate[m,],origin="1970-01-01"),')`s ... Y(',as.Date(est.lag.ydate[m,dim(est.lag.ydate)[2]],origin="1970-01-01"),')`s', sep=""))  
      }
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
  output = list(est.y = est.y,est.ydate = est.ydate, est.x = est.x, est.xdate = est.xdate,
                est.lag.y = est.lag.y, est.lag.ydate = est.lag.ydate,
                out.y = out.y, out.ydate = out.ydate, out.x = out.x, out.xdate = out.xdate,
                out.lag.y = out.lag.y, out.lag.ydate = out.lag.ydate, x.lag = x.lag, y.lag = y.lag,
                min.date = min.date, max.date = max.date)
  return(output)
}
