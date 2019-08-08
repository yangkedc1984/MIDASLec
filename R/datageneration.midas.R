datageneration.midas <- function(rv,options,date=NULL,R=NULL) {
  
  agy <- options$aggrY
  agx <- options$aggrX
  ssize <- length(rv)
  
  
  maxy <- floor((ssize-agx)/agy)
  # LHS
  y <- rowSums(t(matrix(rv[(ssize-maxy*agy+1):ssize],nrow=agy,ncol=maxy,byrow=F))) # <=== sum of realized RV over that horizon agy
  # RHS
  indy <- seq((ssize-maxy*agy),(ssize-agy),by=agy)
  x <- matrix(NA, nrow=length(indy),ncol=agx)
  for (i in 1:length(indy)) {
  x[i,] <- t(rv[seq(indy[i],indy[i]-agx+1,by=-1)])
  }
  Rx <- matrix(NA, nrow=length(indy),ncol=agx)
  if(!is.null(R)){
  for (i in 1:length(indy)) {
    Rx[i,] <- t(R[seq(indy[i],indy[i]-agx+1,by=-1)])
  }
  }
  # dates
  dateselected <- NULL
  if(!is.null(date)){
  begin <- ssize-maxy*agy+1
  indd <- seq(begin,ssize-agy+1,by=agy)
  dateselected <- date[indd] # this gets the first date of the non-overlapping period
  }
  # returns
  Rh <- NULL
  if(!is.null(R)){
  Rh <- rowSums(t(matrix(R[(ssize-maxy*agy+1):ssize],nrow=agy,ncol=maxy,byrow=F)))  # <=== sum of daily log returns over that horizon agy
  }
  return(list(y=y,x=x,dateselected=dateselected,Rh=Rh,Rx=Rx))
}