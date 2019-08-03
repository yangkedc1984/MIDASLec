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
est.start <- "1970-03-01"
est.end <- "1990-03-01"
# mixed.freq.data(data.y, data.ydate, data.x, data.xdate, x.lag, y.lag, horizon, est.start, est.end, disp.flag = TRUE)
data.payems <- mixed.freq.data(rgdp[,2], rgdp[,1], payems[,2], payems[,1], x.lag=9, y.lag=1, horizon=1, est.start=as.Date(est.start), est.end=as.Date(est.end), disp.flag = TRUE)
