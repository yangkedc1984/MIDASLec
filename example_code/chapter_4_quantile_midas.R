# Chapter 4. MIDAS quantile 

# --- preliminaries --- #
rm(list=ls())
require("MIDASLec")
# --- MIDAS quantile regressions --- #

# --- load data --- #
data("example4")

# --- options --- #
optionsmidas <- list(aggrY = 22, aggrX = 22) # forecast horizon, no. of lags of daily returns

# --- compute log returns --- #
snp500[-1, 2] <- log(snp500[-1, 2]/snp500[-dim(snp500)[1], 2])


# --- estimate MIDAS quantile regression with beta (restricted) and plot quantiles --- #

est.midas.0.25 <- midas.optimization.rq(snp500[,2],snp500[,1],optionsmidas,"betaconstr",q.level=0.25,nInitialCond=5,is.plot=TRUE)

est.midas.0.05 <- midas.optimization.rq(snp500[,2],snp500[,1],optionsmidas,"betaconstr",q.level=0.05,nInitialCond=5,is.plot=TRUE)

# --- comppute rolling window unconditional quantile --- #
y <- est.midas.0.05$y
win.size <- 40
un.quant <- matrix(NA,nrow=length(y),ncol=1)
for (j in win.size:length(y)){
  un.quant[j] <- quantile(y[(j-win.size+1):j], probs = 0.05)
}
lines(est.midas.0.05$date,un.quant,type='l',col="red")
# --- estimate CAViAR on monthly data --- #



