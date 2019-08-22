# Chapter 4. MIDAS quantile 

# --- preliminaries --- #
rm(list=ls())
require("MIDASLec")
# --- MIDAS quantile regressions --- #

# --- load data --- #
data("example4")
set.seed(123)

# --- options --- #
optionsmidas <- list(aggrY = 5, aggrX = 22) # forecast horizon, no. of lags of daily returns

# --- compute log returns --- #
snp500[-1, 2] <- log(snp500[-1, 2]/snp500[-dim(snp500)[1], 2])
snp500 <- snp500[-1,]

# --- estimate MIDAS quantile regression with beta (restricted) and plot quantiles --- #
est.midas.0.25 <- midas.optimization.rq(snp500[,2],snp500[,1],optionsmidas,"betaconstr",q.level=0.25,nInitialCond=10)
plot(est.midas.0.25$date,est.midas.0.25$y,type='l',main=paste0("MIDAS quantile. Quantile level: 0.25"), xlab='Weeks',ylab='')
lines(est.midas.0.25$date,est.midas.0.25$cond.quant.MIDAS,type='l',col="red")

est.midas.0.05 <- midas.optimization.rq(snp500[,2],snp500[,1],optionsmidas,"betaconstr",q.level=0.05,nInitialCond=10)
plot(est.midas.0.05$date,est.midas.0.05$y,type='l',main=paste0("MIDAS quantile. Quantile level: 0.05"), xlab='Weeks',ylab='')
lines(est.midas.0.05$date,est.midas.0.05$cond.quant.MIDAS,type='l',col="red")

# --- compute rolling window unconditional quantile --- #
y <- est.midas.0.05$y

win.size <- 100
un.quant <- matrix(NA,nrow=length(y),ncol=1)
for (j in win.size:length(y)){
  un.quant[j] <- quantile(y[(j-win.size+1):j], probs = 0.05)
}
lines(est.midas.0.05$date,un.quant,type='l',col="blue")

# --- estimate CAViAR on weekly data --- #
empiricalQuantile <- un.quant[101] #take unconditional quantile to initialize VaR loop.
date <- est.midas.0.05$date

est.caviar.0.05 <- caviar.optimization(y,q.level=0.05,empiricalQuantile,nInitialCond=10,nInitVec=1000)
plot(date,y,type='l',main=paste0("CAViaR. Quantile level: 0.05"), xlab='Weeks',ylab='')
lines(date,est.caviar.0.05$cond.quant.CAViaR,type='l',col="red")

# --- plot quantiles for MIDAS and CAViaR to compare --- #
plot(date,est.midas.0.05$cond.quant.MIDAS,type='l',main="MIDAS vs CAViaR at 0.05 level", xlab='Weeks',ylab='')
lines(date,est.caviar.0.05$cond.quant.CAViaR,type='l',col="red")

# --- compute conditional skewness at 0.95 level for MIDAS and CAViaR, plot --- #
cond.skewness <- function(cond.up,cond.down,cond.med,q.level){
  num <- cond.up + cond.down - 2 * cond.med
  den <- cond.up - cond.down
  # Kornish-Fisher constant:
  const <- 6/qnorm(q.level)
  c.skew <- num/den*const
  return(c.skew)
}

# --- additionally compute CAViaR for 0.5 and 0.95 levels --- #
est.caviar.0.5 <- caviar.optimization(y,q.level=0.5,quantile(y[1:100],0.5),nInitialCond=10,nInitVec=1000)
est.caviar.0.95 <- caviar.optimization(y,q.level=0.95,quantile(y[1:100],0.95),nInitialCond=10,nInitVec=1000)

# --- compute skewness --- #
skewn.caviar <- cond.skewness(est.caviar.0.95$cond.quant.CAViaR,est.caviar.0.05$cond.quant.CAViaR,est.caviar.0.5$cond.quant.CAViaR,0.95)

# --- additionally compute MIDAS for 0.5 and 0.95 levels --- #
est.midas.0.5  <- midas.optimization.rq(snp500[,2],snp500[,1],optionsmidas,"betaconstr",q.level=0.5,nInitialCond=10)
est.midas.0.95 <- midas.optimization.rq(snp500[,2],snp500[,1],optionsmidas,"betaconstr",q.level=0.95,nInitialCond=10)

skewn.midas <- cond.skewness(est.midas.0.95$cond.quant.MIDAS,est.midas.0.05$cond.quant.MIDAS,est.midas.0.5$cond.quant.MIDAS,0.95)


# --- plot --- #
# --- remove initial 10% due to initialization --- #
idx <- (round(length(skewn.caviar)*0.1)+1):length(skewn.caviar)

plot(date[idx],skewn.midas[idx],type='l',main=paste0("95% conditional skewness estimates"), xlab='Weeks',ylab='', ylim = c(-2, 1))
lines(date[idx],skewn.caviar[idx],type='l',col="red")
