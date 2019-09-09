end.period.date <-
function(date, interval) {
  date.lt <- as.POSIXlt(date) 
  switch(interval, 
         y = {
           date.lt$mon = 11
           date.lt$mday=31
           date=as.Date(date.lt)
         },
         q = {
           date.lt$mon = (date.lt$mon %/% 3 +1)*3 %% 12 
           date.lt$mday = 1
           date.lt$year = date.lt$year + as.integer(date.lt$mon==0)
           date=as.Date(date.lt)-1
         },
         m = {
           date.lt$mon = (date.lt$mon+1) %% 12
           date.lt$mday = 1
           date.lt$year = date.lt$year + as.integer(date.lt$mon==0)
           date=as.Date(date.lt)-1
         },
         d = {
           date = as.Date(date)
         },
         date = as.Date(date$lt)
  )
  return(date)
}
