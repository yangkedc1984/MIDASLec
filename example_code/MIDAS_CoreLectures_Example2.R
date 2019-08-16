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
QLIKE <- RMSFE <- array(NA,c(2,2,5))
for (countf in 1:5) {
  Rd <- AllclosingRet.2018[[countf]][[1]]
  dates <- Alldates.2018[[countf]]
  RV <- AllRV.2018[[countf]][[1]]
  t <- length(dates) # no of time-series obs
  
  
  counth <- 1
  #    loop across horizons
  for (hh in 1:length(Horizons)){
    hhh <- Horizons[hh]
    
    cat("Asset: " , paste0(unlist(Names[[countf]])), "; Horizon: ",  paste0(hhh), '\n')
    optionsmidas$aggrY <- hhh
    
    tmp  <- datageneration.midas(RV,optionsmidas,dates,Rd)
    RVh <- tmp$y
    dates.hhh <- tmp$dateselected
    
    #     MIDAS on RV with beta weights
    beta.est <-  midas.optimization(RV,optionsmidas,"beta")
    # [theta_MIDAS_beta, ~,s2_MIDAS_beta]
    theta.MIDAS.beta <- beta.est$coeff.MIDAS
    s2.MIDAS.beta <- beta.est$cond.val.MIDAS
    #     MIDAS on RV with exp weights
    exp.est <-  midas.optimization(RV,optionsmidas,"exp")
    # [theta_MIDAS_beta, ~,s2_MIDAS_beta]
    theta.MIDAS.exp <- exp.est$coeff.MIDAS
    s2.MIDAS.exp <- exp.est$cond.val.MIDAS
   
    # collect forecasts
    forecasts = cbind(s2.MIDAS.beta,s2.MIDAS.exp)
    
    # evaluate forecasts performance
    qlike <- errors2 <- matrix(NA, nrow=dim(forecasts)[1],ncol=dim(forecasts)[2])
    
    for (fff in 1:2){
        errors2[,fff] <- (RVh - forecasts[,fff])^2
        qlike[,fff] <- log(forecasts[,fff])+RVh/forecasts[,fff]
        RMSFE[counth,fff,countf] <- sqrt(mean(errors2[,fff]))
        QLIKE[counth,fff,countf] <- mean(qlike[,fff])
    }
    
    counth <- counth + 1
  } # end of loop through horizons
  
  
  
} # end of loop through assets