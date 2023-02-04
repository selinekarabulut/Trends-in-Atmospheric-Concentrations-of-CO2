
### Trends in Atmospheric Concentrations of CO2 ###

###load libraries
library(ggfortify)
library(tseries)
library(TSA)
library(GeneCycle)
library(astsa)
library(forecast)

### Data


#get the data
c = read.table("co2.csv", sep=",", header=TRUE)
summary(c)


#plot
co = ts(c[,2])
ts.plot(co, ylab="CO2(PPM)")


#partition data
co2training <- c[1:228, ]
co2clean <- c[1:240, ]
co2full <- c[1:264, ]
#plot data
ts.plot(co2training$co2,ylab="CO2(PPM)")
nt=length(co2training$co2)
fit <- lm(co2training$co2 ~ as.numeric(1:nt)); abline(fit, col="red")
abline(h=mean(co2training$co2), col="blue")


#hist and acf of training data
par(mfrow=c(1,2))
hist(co2training$co2, main="",col="light blue", xlab="")
acf(co2training$co2,lag.max=114,main="")


### Transformations and Stationarity

## Box-Cox
#transformations
# To choose parameter λ of the Box-Cox transformation for datset # Box-Cox transformation:
t = 1:length(co2training$co2)
fit = lm(co2training$co2 ~ t)
bcTransform = boxcox(co2training$co2 ~ t, plotit=TRUE)

lambda=bcTransform$x[which(bcTransform$y == max(bcTransform$y))]
lambda


## Log
par(mfrow=c(1,2))
co2training.bc = (1/lambda)*(co2training$co2^lambda-1)
co2training.log <- log(co2training$co2)
plot.ts(co2training.bc)
plot.ts(co2training.log)


par(mfrow=c(1,3))
hist(co2training$co2, col="light blue", xlab="", main="histogram; CO2 data")
hist(co2training.log, col="light blue", xlab="", main="histogram; ln(U_t)")
hist(co2training.bc, col="light blue", xlab="", main="histogram; bc(U_t)")
#log transform gave a more symmetric histogram and more even variance

#Decomposition of ln(U_t)
library(ggplot2)
library(ggfortify)
y <- ts(as.ts(co2training.log), frequency = 12)
decomp <- decompose(y)
plot(decomp)

#check variance
var(co2training.bc)
var(co2training.log)



## Differencing
par(mfrow=c(1,2))
#differencing at lag 12
co2training.log_12 <- diff(co2training.log, lag=12)
plot.ts(co2training.log_12 , main="log(U_t) differenced at lag 12")
fit <- lm(co2training.log_12~ as.numeric(1:length(co2training.log_12))); abline(fit, col="red")
mean(co2training.log_12)
abline(h=mean(co2training.log_12), col="blue")
#seasonality no longer apparent
# trend is still there
#variance got smaller
#differencing at lag 12 and 1
co2training.stat <- diff(co2training.log_12, lag=1)
plot.ts(co2training.stat, main="log(U_t) differenced at lag 12 and lag 1")
fit <- lm(co2training.stat ~ as.numeric(1:length(co2training.stat))); abline(fit, col="red")
mean(co2training.stat)
abline(h=mean(co2training.stat), col="blue")
#no trend, no seasonality


#check variance
var(co2training.log_12)
var(co2training.stat)

#histogram
hist(co2training.stat, density=20,breaks=20, col="blue", xlab="", prob=TRUE, main="histogram; ln(U_t) differenced at lags 12 & 1")
m<-mean(co2training.stat)
std<- sqrt(var(co2training.stat))
curve( dnorm(x,m,std), add=TRUE )


### Preliminary model identification
#check acf of the transformed and differenced 
acf(co2training.log, lag.max=60, main="ACF of the ln(U_t)")
acf(co2training.log_12, lag.max=60, main="ACF of the ln(U_t), differenced at lag 12")

acf(co2training.stat, lag.max=60, main="")
# ACF outside confidence intervals: Lags 1, maybe 11, 12

#pacf
pacf(co2training.stat, lag.max=60, main="")
# PACF outside confidence intervals: Lags 1, 11, 12, 13, 24, 25,36

# closer inspection
acf(co2training.stat, lag.max=24, main="")



### Model fitting
#List of candidate models to try:s=12, d=1 D=1 p=0,1  P=1 q=1 Q=1

#candidate models
model1 = arima(co2training.log, order=c(0,1,1), seasonal = list(order = c(1,1,1), period = 12), method="ML")
model2=arima(co2training.log, order=c(1,1,1), seasonal = list(order = c(1,1,1), period = 12), method="ML")

# closer inspection
pacf(co2training.stat, lag.max=24, main="")

#check AIC
model1
model2

#revised models
model3 = arima(co2training.log, order=c(0,1,1), seasonal = list(order = c(0,1,1), period = 12), method="ML")
model4=arima(co2training.log, order=c(1,1,1), seasonal = list(order = c(0,1,1), period = 12), method="ML")

#check AIC
model3
model4



### Diagnostic 
#Model 4
par(mfrow=c(1,3))
#To check invertibility of MA part of model 4:
source("plot.roots.R")
plot.roots(NULL,polyroot(c(1, -0.6986)), main="(Model 4) roots of ma part, nonseasonal ")
source("plot.roots.R")
plot.roots(NULL,polyroot(c(1, -0.9169)), main="(Model 4) roots of ma part, seasonal ")
#To check stationarity of AR part of model 4:
source("plot.roots.R")
plot.roots(NULL,polyroot(c(1, -0.2766)), main="(Model 4) roots of ar part, nonseasonal")

#roots of AR
polyroot(c(1, -0.2766))

#roots of MA
polyroot(c(1, -0.6986))

#roots of seasonal MA
polyroot(c(1, -0.9169))

#Roots of AR, MA and seasonal MA part lie outside the unit circle. 
#Hence, model 4 is both stationary and invertible. I will follow with series of diagnostics and their interpretation below.

reslog4 <- residuals(model4)
#plot of residuals
plot.ts(reslog4)
fitres4 <- lm(reslog4 ~ as.numeric(1:length(reslog4))); abline(fitres4, col="red")
abline(h=mean(reslog4), col="blue")

par(mfrow=c(1,2))
#histogram of residuals
hist(reslog4,density=20,breaks=40, col="blue", xlab="", prob=TRUE,main="",cex=1)
m4 <- mean(reslog4)
std4 <- sqrt(var(reslog4))
curve( dnorm(x,m4,std4), add=TRUE )
#qq plot of residuals
qqnorm(reslog4,main= "",cex=1)
qqline(reslog4,col="blue")

par(mfrow=c(1,2))
#acf, pacf of model 4
acf(reslog4, lag.max=60,
    main="",cex=1)
pacf(reslog4, lag.max=60,
     main="",cex=1)

shapiro.test(reslog4)

Box.test(reslog4, lag = 15, type = c("Box-Pierce"), fitdf = 3)

Box.test(reslog4, lag = 15, type = c("Ljung-Box"), fitdf = 3)

Box.test((reslog4)^2, lag = 15, type = c("Ljung-Box"), fitdf = 0)


#Model 3
par(mfrow=c(1,2))
#To check invertibility of MA part of model 3:
source("plot.roots.R")
plot.roots(NULL,polyroot(c(1, -0.4669)), main="(Model 3) roots of ma part, nonseasonal ")
source("plot.roots.R")
plot.roots(NULL,polyroot(c(1, -0.9026)), main="(Model 3) roots of ma part, seasonal ")

reslog3 <- residuals(model3)
#plot of residuals
plot.ts(reslog3)
fitres3 <- lm(reslog3 ~ as.numeric(1:length(reslog3))); abline(fitres3, col="red")
abline(h=mean(reslog3), col="blue")

par(mfrow=c(1,2))
#histogram of residuals
hist(reslog3,density=20,breaks=40, col="blue", xlab="", prob=TRUE,main="",cex=1)
m3 <- mean(reslog3)
std3 <- sqrt(var(reslog3))
curve( dnorm(x,m3,std3), add=TRUE )
#qq plot of residuals
qqnorm(reslog3,main= "",cex=1)
qqline(reslog3,col="blue")

par(mfrow=c(1,2))
#acf, pacf of model 3
acf(reslog3, lag.max=60,
    main="",cex=1)
pacf(reslog3, lag.max=60,
     main="",cex=1)

shapiro.test(reslog3)

Box.test(reslog3, lag = 15, type = c("Box-Pierce"), fitdf = 2)

Box.test(reslog4, lag = 15, type = c("Ljung-Box"), fitdf = 2)

Box.test((reslog4)^2, lag = 15, type = c("Ljung-Box"), fitdf = 0)



### Forecasting
co = ts(c[,2])
co2train1 = co[c(1:228)]
co2clean1 =co[c(1:240)]
co2full1 = co[c(1:264)]

#forecast on log
pred.tr <- predict(model4, n.ahead = 12)
U.tr= pred.tr$pred + 2*pred.tr$se
L.tr= pred.tr$pred - 2*pred.tr$se
#forecast on original
pred.orig=exp(pred.tr$pred)
U=exp(U.tr)
L=exp(L.tr)
#plot-zoomed
ts.plot(co2clean1, xlim = c(200,length(co2train1)+12), ylim = c(400,max(U)), col="red",
        main="",
        ylab="CO2(PPM")
lines(U, col="blue", lty="dashed")
lines(L, col="blue", lty="dashed")
points((length(co2train1)+1):(length(co2train1)+12), pred.orig, col="black")


pred4 <- predict(model4, n.ahead = 36)
U4= pred4$pred + 2*pred4$se
L4= pred4$pred - 2*pred4$se
pred4.orig=exp(pred4$pred)
U4orig=exp(U4)
L4orig=exp(L4)
ts.plot(co2full1, xlim = c(200,length(co2train1)+36), ylim = c(400,max(U4orig)), col="red",lwd=2,
        main="",
        ylab="CO2(PPM)")
lines(U4orig, col="blue", lty="dashed")
lines(L4orig, col="blue", lty="dashed")
points((length(co2train1)+1):(length(co2train1)+36), pred4.orig, col="black")



### Spectral Analysis
library(TSA)
TSA::periodogram(reslog4,main="")

fisher.g.test(reslog4)

cpgram(reslog4,main="")














