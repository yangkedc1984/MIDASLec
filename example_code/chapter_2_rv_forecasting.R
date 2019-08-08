# Chapter 2. RV forecasting

# --- preliminaries --- #
rm(list=ls())
require("MIDASLec")
# --- MIDAS RV regressions --- #

# --- load data --- #
data("example2")

# --- options --- #
optionsmidas <- list(aggrX =  126) # no. of lags of RV
#Horizons = c(5, 10, 22, 44, 66) # forecasting horizons
Horizons = c(5, 10) # forecasting horizons


# --- main loop through assets --- #

for (countf in 1:5) {
  Rd <- AllclosingRet.2018[[countf]]
  dates <- Alldates.2018[[countf]]
  RV <- AllRV_2018{countf,1}
  t <- size(dates,1) # no of time-series obs
  
  
  
  
} # end of loop through assets