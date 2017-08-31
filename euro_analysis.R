library(xts)
library(fBasics)
library(matrixStats)
library(fGarch)

setwd("~/Desktop/Time Series/final")
load("/Users/guanwei/Desktop/Time Series/final/alldata.RData")
statfun<-function(x,func_name){
  func_name <- match.fun(func_name)
  m<-xts(rep(0,nrow(x)),order.by=index(x))
  for(i in 1:nrow(x)){
    idx2<-col(x[i,])[is.na(x[i,])==TRUE]
    p<-x[i,-idx2]
    idx<-col(p)[duplicated(p)==TRUE]
    m[i,]<-func_name(p[,-idx])
  } 
  m
}

assetmed.eur<-statfun(asset_eur,median)
assetsum.eur<-statfun(asset_eur,sum)
assetmean.eur<-statfun(asset_eur,mean)
assetsd.eur<-statfun(asset_eur,sd)

plot.zoo(merge(assetsum.eur,assetmed.eur,assetmean.eur,assetsd.eur))

wammed.eur<-statfun(wam_eur,median)
wamsum.eur<-statfun(wam_eur,sum)
wammean.eur<-statfun(wam_eur,mean)
wamsd.eur<-statfun(wam_eur,sd)

plot.zoo(merge(wamsum.eur,wammed.eur,wammean.eur,wamsd.eur))

# walmed.eur<-statfun(wal_eur,median)
# walsum.eur<-statfun(wal_eur,sum)
# walmean.eur<-statfun(wal_eur,mean)
# walsd.eur<-statfun(wal_eur,sd)
# 
# plot.zoo(merge(walsum.eur,walmed.eur,walmean.eur,walsd.eur))


#arima
plot(assetsum.eur)
dassetsum.eur=diff(log(assetsum.eur))
plot(dassetsum.eur)
dassetsum.eur.ts=ts(dassetsum.eur,start=c(2007,1,5),freq=52)
dassetsum.eur.ts=na.omit(dassetsum.eur.ts)
f1 <- acf(dassetsum.eur.ts)
f1$lag <- f1$lag*52
plot(f1,xlab="Lag(weeks)")
g1 <- pacf(dassetsum.eur.ts)
g1$lag <- g1$lag*52
plot(g1,xlab="Lag(Weeks)")
g1$acf
t.test(dassetsum.eur)

Box.test(dassetsum.eur,lag=13,type = 'Ljung')
m1=ar(dassetsum.eur.ts,method = 'mle')
m1$order
m1=arima(dassetsum.eur.ts,order=c(13,0,0))
m1
tratio=m1$coef/sqrt(diag(m1$var.coef))
tratio
m1=arima(dassetsum.eur.ts,order=c(13,0,0),fixed = c(NA,0,0,NA,0,0,0,0,NA,0,0,0,NA,NA))
m1
mm1<-arima(dassetsum.eur.ts,c(13,0,0),seasonal = list(order= c(0,0,1),period=10))
mm1
tratio=mm1$coef/sqrt(diag(mm1$var.coef))
tratio
mm1<-arima(dassetsum.eur.ts,c(13,0,0),seasonal = list(order= c(0,0,1),period=10),fixed = c(NA,0,0,NA,0,0,0,0,NA,NA,0,NA,NA,NA,0))
mm1
tratio=mm1$coef/sqrt(diag(mm1$var.coef))
tratio
tsdiag(m1,gof=12)
Box.test(m1$residuals,lag = 12,type = 'Ljung')


dwammed.eur=diff(log(wammed.eur))
plot(dwammed.eur)
t.test(dwammed.eur)
dwammed.eur.ts=ts(dwammed.eur,start=c(2007,1,5),freq=52)
dwammed.eur.ts=na.omit(dwammed.eur.ts)
f2 <- acf(dwammed.eur.ts)
f2$lag <- f2$lag*52
plot(f2,xlab="Lag(weeks)")
g2 <- pacf(dwammed.eur.ts)
g2$lag <- g2$lag*52
plot(g2,xlab="Lag(Weeks)")
g2$acf

Box.test(dwammed.eur.ts,lag=12,type = 'Ljung')
m2=ar(dwammed.eur.ts,method = 'mle')
m2$order
m2=arima(dwammed.eur.ts,order=c(8,0,6))
m2
tratio=m2$coef/sqrt(diag(m2$var.coef))
tratio
m2=arima(dwammed.eur.ts,order=c(8,0,6),fixed = c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0,0))
m2
mm2<-arima(dwammed.eur.ts,c(8,0,6),seasonal = list(order= c(0,0,1),period=10))
mm2
tratio=mm2$coef/sqrt(diag(mm2$var.coef))
tratio
mm2<-arima(dwammed.eur.ts,c(8,0,6),seasonal = list(order= c(0,0,1),period=10),fixed = c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0,NA,0))
mm2

tsdiag(m2,gof=12)
Box.test(m2$residuals,lag = 12,type = 'Ljung')

#arch & Garch
v <- dassetsum.eur.ts-mean(dassetsum.eur.ts)
Box.test(v^2,lag=12,type="Ljung-Box")

u <- dwammed.eur.ts-mean(dwammed.eur.ts)
Box.test(u^2,lag=12,type="Ljung-Box")
Box.test(dwammed.eur.ts,lag=12,type="Ljung-Box")
archTest=function(rtn,m=10){
  # Perform Lagrange Multiplier Test for ARCH effect of a time series
  # rtn: time series
  # m: selected AR order
  #
  y=(rtn-mean(rtn))^2
  T=length(rtn)
  atsq=y[(m+1):T]
  x=matrix(0,(T-m),m)
  for (i in 1:m){
    x[,i]=y[(m+1-i):(T-i)]
  }
  md=lm(atsq~x)
  summary(md)
}

archTest(u,12)
#F=94.5 with p-value of 2.2*10^(-16). The test confirms strong ARCH effects in the monthly log returns
archTest(v,12)

m3 <- garchFit(~arma(5,4)+garch(1,1),data=dassetsum.eur.ts,trace=F)
summary(m3)
resi2 <- residuals(m3,standardize=T)
pacf(resi2,lag=20)
pacf(resi2^2,lag=20)
plot(m3)
13
0

m4 <- garchFit(~arma(3,2)+garch(1,1),data=dwammed.eur.ts,trace=F)
summary(m4)
resi2 <- residuals(m4,standardize=T)
pacf(resi2,lag=20)
pacf(resi2^2,lag=20)
plot(m4)
13
0

tdx=(1:448)/52+2007

source("backtest.R")
pm1=backtest(m1,dassetsum.eur.ts,400,1)
pm1fit=dassetsum.eur.ts[401:448]-pm1$error
plot(tdx[401:448],dassetsum.eur.ts[401:448],xlab='year',ylab='asset growth',type='l')
points(tdx[401:448],pm1fit,pch='*')

pm2=backtest(m2,dwammed.eur.ts,400,1)
pm2fit=dwammed.eur.ts[401:448]-pm2$error
plot(tdx[401:448],dwammed.eur.ts[401:448],xlab='year',ylab='wam growth',type='l')
points(tdx[401:448],pm2fit,pch='*')

source("backtestGarch.R")
pm3=backtestGarch(dassetsum.eur.ts,400,1,inc.mean = F,cdist = "sstd")
pm4=backtestGarch(dwammed.eur.ts,400,1,inc.mean = F,cdist = "sstd")



