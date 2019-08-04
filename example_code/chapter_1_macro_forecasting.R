# Chapter 1. Macro forecasting

# --- preliminaries --- #
rm(list=ls())
require("MIDASLec")
# --- ADL-MIDAS(p,q) regressions --- #

# --- load data --- #
data("example1")

# --- check if data is loaded in your RStudio (Environment) --- #
# --- you should see: cfnai, payems, rgdp --- # 

# --- trasform data to growth rates --- #
# --- cfnai is the first PC, no need to transform --- #
# --- payems: --- # 
payems[-1, 2] <- log(payems[-1, 2]/payems[-dim(payems)[1], 2])*100
payems <- payems[-1, ]
#plot(payems[,1],payems[,2],type='l',ylab='Monthly employment growth rate',xlab='Months')
# --- rgdp: --- #
rgdp[-1, 2] <- log(rgdp[-1, 2]/rgdp[-dim(rgdp)[1], 2])*100
rgdp <- rgdp[-1, ]
#plot(rgdp[,1],rgdp[,2],type='l',ylab='Quarterly real GDP growth rate',xlab='Quarters')


# --- construct MIDAS data structures to later use in models --- #
# --- First example is with Employment data: --- #
# --- intial and last date for in-sample estimation --- #
est.start <- "1985-01-01"
est.end <- "2009-01-01"

# mixed.freq.data(data.y, data.ydate, data.x, data.xdate, x.lag, y.lag, horizon, est.start, est.end, disp.flag = TRUE)
data.payems.in <- mixed.freq.data(rgdp[,2], rgdp[,1], payems[,2], payems[,1], x.lag=9, y.lag=1, horizon=3, est.start=as.Date(est.start), est.end=as.Date(est.end), disp.flag = TRUE)


# --- Estimate ADL-MIDAS regression model using employment data --- #
weight <- nealmon
# --- get initial values --- #
set.seed(123)
# startx.all <- get.start.adl.midas(y=data.payems.in$est.y,X=data.payems.in$est.x,z=data.payems.in$est.lag.y,weight=weight,par.num.weight=3,num.evals=10000, num.coef=100)
startx.all <- get.start.adl.midas(y=data.payems.in$est.y,X=data.payems.in$est.x,z=data.payems.in$est.lag.y,weight=weight,par.num.weight=3,num.evals=10000, num.coef=1)

# est <- vals <- NULL
# for (j in 1:dim(startx.all)[1]) { 
#   est[[j]] <- midas_r_plain(y=data.payems.in$est.y,X=data.payems.in$est.x,z=cbind(data.payems.in$est.lag.y,rep(1,times=length(data.payems.in$est.y))),weight=weight,startx=startx.all[j,-((dim(startx.all)[2]-1):dim(startx.all)[2])],startz=startx.all[j,((dim(startx.all)[2]-1):dim(startx.all)[2])],control=list(maxit=500))
#   vals[j] <- est[[j]]$opt$value[1]
# }
#est.midas <- est[[which(min(vals)==vals)]]
est.midas <- midas_r_plain(y=data.payems.in$est.y,X=data.payems.in$est.x,z=cbind(data.payems.in$est.lag.y,rep(1,times=length(data.payems.in$est.y))),weight=weight,startx=startx.all[-c(4,5)],startz=startx.all[c(4,5)],control=list(maxit=500))

# --- Estimate ADL regression model using time-averaged employment data --- # 
est.ta <- lm(data.payems.in$est.y~data.payems.in$est.lag.y+rowMeans(data.payems.in$est.x[,1:3])+rowMeans(data.payems.in$est.x[,4:6])+rowMeans(data.payems.in$est.x[,7:9]))


# --- compare IS --- #
sqrt(mean(est.midas$residuals^2))
sqrt(mean(est.ta$residuals^2))

# --- compare OOS --- #
midas.oos <- forecast.adl(est.midas,weight=weight,par.num.weight=3,data.payems.in,is.intercept=TRUE)
ta.oos <- forecast.ta.example(est.ta,data.payems.in)

midas.oos$rmse
ta.oos$rmse


# --- plot lag polynomials --- #
par(mfrow=c(1,2)) 
plot(nealmon(est.midas$coefficients[c(1,2,3)],d=9),type='l',
      xlab='Lag',ylab='Coefficient',
      main='Normalized exponential Almon lag polynomial')

plot(c(rep(as.numeric(est.ta$coefficients[3]),times=3),rep(as.numeric(est.ta$coefficients[4]),times=3),rep(as.numeric(est.ta$coefficients[5]),times=3)),type='l',
      xlab='Lag',ylab='Coefficient',
      main='Time-averaged data coefficients')


# --- compare fixed, rolling, expanding windows for predictions --- #
# --- fixed --- #
midas.obj.fixed <- midas.adl(data.y=rgdp[,2],data.ydate=rgdp[,1],
                             data.x=payems[,2],data.xdate=payems[,1], 
                             est.start=as.Date(est.start),est.end=as.Date(est.end),
                             horizon=3,x.lag=9,y.lag=1,polynomial="nealmon",method="fixed",disp.flag=TRUE,
                             num.evals=10000,num.coef=1)

ta.obj.fixed <- midas.adl(data.y=rgdp[,2],data.ydate=rgdp[,1],
                          data.x=payems[,2],data.xdate=payems[,1], 
                          est.start=as.Date(est.start),est.end=as.Date(est.end),
                          horizon=3,x.lag=9,y.lag=1,polynomial="timeaverage",method="fixed",disp.flag=TRUE)
# --- rolling --- #
midas.obj.rolling <- midas.adl(data.y=rgdp[,2],data.ydate=rgdp[,1],
                               data.x=payems[,2],data.xdate=payems[,1], 
                               est.start=as.Date(est.start),est.end=as.Date(est.end),
                               horizon=3,x.lag=9,y.lag=1,polynomial="nealmon",method="rolling",disp.flag=TRUE,
                               num.evals=10000,num.coef=1)

ta.obj.rolling <- midas.adl(data.y=rgdp[,2],data.ydate=rgdp[,1],
                            data.x=payems[,2],data.xdate=payems[,1], 
                            est.start=as.Date(est.start),est.end=as.Date(est.end),
                            horizon=3,x.lag=9,y.lag=1,polynomial="timeaverage",method="rolling",disp.flag=TRUE)
# --- expanding --- #
midas.obj.expand <- midas.adl(data.y=rgdp[,2],data.ydate=rgdp[,1],
                              data.x=payems[,2],data.xdate=payems[,1], 
                              est.start=as.Date(est.start),est.end=as.Date(est.end),
                              horizon=3,x.lag=9,y.lag=1,polynomial="nealmon",method="expand",disp.flag=TRUE,
                              num.evals=10000,num.coef=1)

ta.obj.expand <- midas.adl(data.y=rgdp[,2],data.ydate=rgdp[,1],
                           data.x=payems[,2],data.xdate=payems[,1], 
                           est.start=as.Date(est.start),est.end=as.Date(est.end),
                           horizon=3,x.lag=9,y.lag=1,polynomial="timeaverage",method="rolling",disp.flag=TRUE)

# --- Second example is with CFNAI data: --- #
# --- intial and last date for in-sample estimation --- #
est.start <- "1987-01-01" 
est.end <- "2011-12-01"
data.cfani.in <- mixed.freq.data(rgdp[,2], rgdp[,1], cfnai[,2], cfnai[,1], x.lag=12, y.lag=1, horizon=3, est.start=as.Date(est.start), est.end=as.Date(est.end), disp.flag = TRUE)

# --- Estimate ADL-MIDAS regression model using CFNAI data --- #
weight <- nbeta
# --- get initial values --- #
set.seed(123)
# startx.all <- get.start.adl.midas(y=data.payems.in$est.y,X=data.payems.in$est.x,z=data.payems.in$est.lag.y,weight=weight,par.num.weight=3,num.evals=10000, num.coef=100)
startx.all <- get.start.adl.midas(y=data.cfani.in$est.y,X=data.cfani.in$est.x,z=data.cfani.in$est.lag.y,weight=weight,par.num.weight=3,num.evals=10000, num.coef=1)


est.midas <- midas_r_plain(y=data.cfani.in$est.y,X=data.cfani.in$est.x,z=cbind(data.cfani.in$est.lag.y,rep(1,times=length(data.cfani.in$est.y))),weight=weight,startx=startx.all[-c(4,5)],startz=startx.all[c(4,5)],control=list(maxit=500))

# --- Estimate ADL regression model using U-MIDAS scheme --- # 
est.umidas <- lm(data.cfani.in$est.y~data.cfani.in$est.lag.y+data.cfani.in$est.x)


# --- compare IS --- #
sqrt(mean(est.midas$residuals^2))
sqrt(mean(est.umidas$residuals^2))

# --- compare OOS --- #
midas.oos <- forecast.adl(est.midas,weight=weight,par.num.weight=3,data.cfani.in,is.intercept=TRUE)
umidas.oos <- forecast.umidas(est.umidas,data.cfani.in)

midas.oos$rmse
umidas.oos$rmse

# --- plot lag polynomials --- #
par(mfrow=c(1,2)) 
plot(nbeta(est.midas$coefficients[c(1,2,3)],d=12),type='l',
     xlab='Lag',ylab='Coefficient',
     main='Normalized Beta lag polynomial')

plot(coef(est.umidas),type='l',
     xlab='Lag',ylab='Coefficient',
     main='U-MIDAS lag polynomial')


