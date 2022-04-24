# CHANGE POINT DETECTION

################################################################################
## This task is an attempt to check “whether the pattern of the Wells Fargo 
## stock follows a Geometric Brownian Motion or not” due to detecting changes. 
## Answering this question needs to check for either independency or 
## normality of the log returns in the presence of Black & Scholes assumption.
################################################################################

## install.packages("pacman")
## Load required packages at once
library(pacman)
pacman::p_load(quantmod,yuima,sde,fOptions)

library(quantmod)
library(yuima)
library(sde)
library(fOptions)

## Check (in)dependency of the log returns through correlation analysis


## Load Data from Yahoo Finance and make data preparation
getSymbols("WFC", from="2014-01-01")
s<-WFC$WFC.Adjusted
x<-diff(log(s))
x<-na.omit(x)

## Plot autocorrelation of log returns 
acf(x,main="ACF plot for Log Returns", col="red")

## Check Normality by p-value through Shapiro-Test
shapiro.test(as.numeric(x))

## Assess abnormality for the log returns with Graphical Q-Q plot
par(mfrow=c(2,1))
y<-rnorm(10000)
qqnorm(y)
qqline(y, col="red")
qqnorm(x)
qqline(x, col="red")

## Show the discrepancy between the theoretical Gaussian distribution and
## the distribution followed by the log returns of the Wells Fargo stock
hist(x, breaks=30, freq=FALSE, xlab = "log return", main="DENSITY vs. NORMAL DISTRIBUTION")
f<-function(y) dnorm(y, mean=mean(x), sd=sd(x))
curve(f, min(x),max(x), col="blue", lwd=3, add=TRUE)
lines( density(x), col="red",lwd=3 )

## Estimate the percentage drift and the percentage volatility
## with assumption of Geometric Brownian motion (GBM)
Dt <- 1/252 # Data Frequency
alpha.hat <- mean(x)/Dt
variance.hat <- var(x)/Dt
variance.hat
sigma.hat <- sqrt(variance.hat)
sigma.hat
mu.hat<-alpha.hat+0.5* variance.hat
mu.hat

## Depict possible pattern of whatever stock following a GBM behavior
mod <- setModel(drift="0.105*x", diff="0.15629 *x", solve.variable = c("x"))
samp <- setSampling(T=10, n=100000)
x <- simulate(mod, sampling=samp,xinit=5)
plot(x,col="blue",main="GBM PATTERN SIMUlATION", xlab="Time", ylab="Values")

## Perform model selection process when stock price does not follow 
## a Geometric Brownian Motion,:
## To do this, use continuous models using 
## Akaike Information Criterion (AIC) statistics
## e.g. Merton Model; Vasicek Model;
## Cox, Ingersoll & Ross Models; Dothan Model; Geometric Brownian Motion Model;
## Brennan & Schwartz Model; Constant Elasticity Variance Model; and CKLS Model.

## Load and Prepare Data
getSymbols("WFC", from="2014-01-01")
s<-WFC$WFC.Adjusted
X<-diff(log(s))
X<-na.omit(X)
m <- mean(s)
sigma<- as.numeric(sqrt(var(X)*252)) # n=252, the number of working days in a year
M <- mean(X)*252

# Create facilitator matrices containing the parameters of the models
drift <- c("mu*x", "alpha", "alpha*(mu-x)", "0", "alpha*(mu-x)", 
           "alpha*(mu-x)", "0", "alpha*x", "alpha*(mu-x)")
diffusion <- c("sigma*x", "sigma", "sigma*sqrt(x)", "sigma*x", "sigma", 
               "sigma*x", "sigma*x^(3/2)", "sigma*x^gamma", "sigma*x^gamma")
AIC <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA)
MyMat <- rbind(drift,diffusion,AIC)
colnames(MyMat) <- c("GBM", "Merton", "CIR", "Dothan", "Vasicek", "BrenSch", "CIR2", "CEV", "CKLS")
row.names(MyMat) <- c("drift", "diffusion", "AIC")

qlme_GBM <- c(NA, M, sigma,NA); qlme_Merton <- c(1, NA, sigma, NA); qlme_CIR <- c(1, m, sigma, NA)
qlme_Dothan <- c(NA, NA, sigma, NA); qlme_Vasicek <- c(1, m, sigma, NA); qlme_BrenSch <- c(1, m, sigma, NA)
qlme_CIR2 <- c(NA, NA, sigma, NA); qlme_CEV <- c(1, NA, sigma, 0.5); qlme_CKLS <- c(1, m, sigma, 0.5)
qlme_matrix <- rbind(qlme_GBM, qlme_Merton, qlme_CIR, qlme_Dothan, qlme_Vasicek, qlme_BrenSch, qlme_CIR2, qlme_CEV, qlme_CKLS)
colnames(qlme_matrix) <- c("alpha", "mu", "sigma", "gamma")
row.names(qlme_matrix) <- c("GBM", "Merton", "CIR", "Dothan", "Vasicek", "BrenSch", "CIR2", "CEV", "CKLS")

# Find the preferred model with the lowest AIC value.

for (i in 1:dim(MyMat)[2]) {
  SMoD <- setModel(drift=MyMat[1,i], diff=MyMat[2,i], solve.variable = c("x"))
  yuima.SMoD <- setYuima(model=SMoD, data=setData(s, delta=1/252))
  if (i == 1){ 
    fit.SMoD <- qmle(yuima.SMoD, start=list(mu=qlme_matrix[i,2], sigma=qlme_matrix[i,3]))
  } else if(i==2){
    fit.SMoD <- qmle(yuima.SMoD, start=list(alpha=qlme_matrix[i,1], sigma=qlme_matrix[i,3]))
  } else if(i==3){
    fit.SMoD <- qmle(yuima.SMoD, start=list(alpha=qlme_matrix[i,1], mu=qlme_matrix[i,2], sigma=qlme_matrix[i,3]))
  } else if(i==4){
    fit.SMoD <- qmle(yuima.SMoD, start=list(sigma=qlme_matrix[i,3]))
  }else if(i==5 || i==6){
    fit.SMoD <- qmle(yuima.SMoD, start=list(alpha=qlme_matrix[i,1], mu=qlme_matrix[i,2], sigma=qlme_matrix[i,3]))
  }else if(i==7){
    fit.SMoD <- qmle(yuima.SMoD, start=list(sigma=qlme_matrix[i,3]))
  }else if(i==8){
    fit.SMoD <- qmle(yuima.SMoD, start=list(alpha=qlme_matrix[i,1], sigma=qlme_matrix[i,3], gamma=qlme_matrix[i,4]))
  }else {
    fit.SMoD <- qmle(yuima.SMoD, start=list(alpha=qlme_matrix[i,1], mu=qlme_matrix[i,2], sigma=qlme_matrix[i,3], gamma=qlme_matrix[i,4]))
  }
  MyMat[3,i] <- AIC(fit.SMoD)
}

which.min(MyMat[3,])
MyMat[3,which.min(MyMat[3,])]


## Utilizing Change Point Analysis try to detect 
## volatility clustering in stock market

## Load and order data
getSymbols("WFC", from="2014-01-01", to="2015-08-25")
S <- zoo(WFC$WFC.Close)

## Find and show the volatility change-point
cpoint(S)
plot(S, xlab="TIME", ylab="PRICE", main="STOCK PRICE PLOT AND CHANGE POINT ANALYSIS")
abline(v=cpoint(S), col="red",lwd=2)


## Compute the differences among estimated values 
## by using the selected model with respect to the speed of 
## mean reversion, long range mean and volatility.

## Load and Prepare data
getSymbols("WFC", from="2014-01-01", to="2015-08-25")
s<-zoo(WFC$WFC.Adjusted)
X<-diff(log(s))
X<-na.omit(X)
m <- mean(s)
sigma<- as.numeric(sqrt(var(X)*252))
getSymbols("WFC", from="2014-01-01", to="2014-10-05")
s1<-zoo(WFC$WFC.Adjusted)
getSymbols("WFC", from="2014-10-06", to="2015-08-25")
s2<-zoo(WFC$WFC.Adjusted)

## Calculate the quasi-likelihood and estimate of the parameters of 
## the stochastic differential equation by 
## the maximum likelihood method or 
## least squares estimator of the drift parameter.

## Build the stochastic differential equations
brensch <- setModel(drift="alpha*(mu-x)", diffusion="sigma*x", solve.variable = c("x"))

## Construct an object of by combining "model", "data",
## "sampling", "characteristic" and "functional"slots.
yuima <- setYuima(model=brensch, data=setData(s, delta=1/252))
yuima1 <- setYuima(model=brensch, data=setData(s1, delta=1/252))
yuima2 <- setYuima(model=brensch, data=setData(s2, delta=1/252))
start <- list(alpha=1, mu=m, sigma=sigma)

## Calculate the quasi-likelihood and ML estimator
fit <- qmle(yuima, start=start)
summary(fit)

