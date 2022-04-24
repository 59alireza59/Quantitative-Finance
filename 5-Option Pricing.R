# OPTION PRICING

## This task represents a framework to compute option price 
## using Black and Scholes formula and Monte Carlo simulation 
## in which the first step requirement is to apply Changing Point Analysis
## over a rich financial time series dataset from Yahoo Finance (Johnson &
## Johnson stock prices). The changing point analysis helps to find when 
## volatility changes in data with respect to the length of time considered.

## install.packages("pacman")
## Load required packages at once
library(pacman)
pacman::p_load(tseries,sde,quantmod,fOptions)

library(quantmod)
library(tseries)
library(sde)
library(fOptions)

## Changing Point Analysis
CP1 <- get.hist.quote("JNJ",start="2014-01-01", end = "2015-05-05")
CP1 <- CP1$Close
cpoint(CP1)

CP2 <- get.hist.quote("JNJ",start="2014-09-30", end = "2015-05-05")
CP2 <- CP2$Close
cpoint(CP2)

CP3 <- get.hist.quote("JNJ",start="2014-10-17", end = "2015-05-05")
CP3 <- CP3$Close
cpoint(CP3)

CP4 <- get.hist.quote("JNJ",start="2014-12-09", end = "2015-05-05")
CP4 <- CP4$Close
cpoint(CP4)

plot(CP1)
abline(v=cpoint(CP1)$tau0,lty=3,col="green")
abline(v=cpoint(CP2)$tau0,lty=3,col="red")
abline(v=cpoint(CP3)$tau0,lty=3,col="blue")
abline(v=cpoint(CP4)$tau0,lty=3,col="black")

## Pricing with Black and Scholes formula / GBM

## Load and Prepare Data
getSymbols("JNJ",from="2014-09-30",to="2015-05-05")
S <- JNJ$JNJ.Close
plot(S,main="JNJ.Adjusted")
X <- diff(log(S))
X <- na.omit(X)
plot(X,main="JNJ.Return")

## Black and Scholes Formulation
Delta <- 1/252
alpha.hat <- mean(X,na.rm=TRUE)/Delta
sigma.hat <- sqrt(var(X,na.rm=TRUE)/Delta)
mu.hat <- alpha.hat+0.5*sigma.hat^2
mu.hat
sigma.hat
alpha.hat
S0 <- 99.50
r <- 0.03
sigma <- 0.1739121
T <- 19/252

## Calculate the option price using Call and Put Options
## under assumption of Black and Scholes theory
K <- c(98, 102, 105)
BS_Call <-matrix(0,nrow=3,ncol=1)
colnames(BS_Call) <- c("BS_Call"); row.names(BS_Call) <- c(98, 102, 105)
BS_Put <- matrix(0,nrow=3,ncol=1)
colnames(BS_Put) <- c("BS_Put"); row.names(BS_Put) <- c(98, 102, 105)
for (i in 1:length(K)) {
  BS_Call[i,1] <- GBSOption(TypeFlag = "c", S = S0, X = K[i], Time = T, r = r, b = r, sigma = sigma)@price
  BS_Put[i,1] <- GBSOption(TypeFlag = "p", S = S0, X = K[i], Time = T, r = r, b = r, sigma = sigma)@price
}

## Look at Implied / Smiles volatility

## Get Call option
sigma.hat <- 0.1739121
KC <- c(2.47,0.27,0.07)
K <- c(98,102,105)
S0 <- 99.5
nP <- length(KC)
T <- 19/252
r <- 0.03

## Get Smiles volatility with Call options
smile <- sapply(1:nP, function(i) GBSVolatility(KC[i], "c", S = S0, X = K[i], Time = T, r = r, b = r))
vals <- c(smile, sigma.hat)
plot(K, smile, type = "l", ylim = c(min(vals, na.rm = TRUE), max(vals, na.rm = TRUE)), main 
     ="Volatility smile: Call")
abline(v = S0, lty = 3, col = "blue")
abline(h = sigma.hat, lty = 3, col = "red")
axis(2, sigma.hat, expression(hat(sigma)), col = "red")


# Get Put option
sigma.hat <- 0.1739121
KP <- c(1.09,2.59,5.2)
K <- c(98,102,105)
S0 <- 99.5
nP <- length(KP)
T <- 19/252
r <- 0.03

## Get Smiles volatility with Put options
smile <- sapply(1:nP, function(i) GBSVolatility(KP[i], "p", S = S0, X = K[i], Time = T, r = r, b = r))
vals <- c(smile, sigma.hat)
plot(K, smile, type = "l", ylim = c(min(vals, na.rm = TRUE), max(vals, na.rm = TRUE)), main = "Volatily smile: Put")
abline(v = S0, lty = 3, col = "blue")
abline(h = sigma.hat, lty = 3, col = "red")
axis(2, sigma.hat, expression(hat(sigma)), col = "red")

## Monte Carlo Method Option Pricing

## Make Monte Carlo Price
S0 <- 99.50
r <- 0.03
sigma <- 0.1739121
T <- 19/252

MCPrice<-function(x=1,t=0,T=1,r=1,sigma=1,M=1000,f){
  h<-function(m){
    u<-rnorm(m/2)
    tmp<-c(x*exp((r-0.5*sigma^2)*(T-t)+sigma*sqrt(T-t)*u),
           x*exp((r-0.5*sigma^2)*(T-t)+sigma*sqrt(T-t)*(-u)))
    mean(sapply(tmp,function(xx) f(xx)))
  }
  p<-h(M)
  p*exp(-r*(T-t))
}

## Look at Call and Put options using Make Monte Carlo simulation
K <- c(98, 102, 105)
MC_Call <-matrix(0,nrow=3,ncol=3)
colnames(MC_Call) <- c("MC_Call1","MC_Call2","MC_Call3")
row.names(MC_Call) <- c(98, 102, 105)
MC_Put <- matrix(0,nrow=3,ncol=3)
colnames(MC_Put) <- c("MC_Put1","MC_Put2","MC_Put3")
row.names(MC_Put) <- c(98, 102, 105)
N_Simul <- c(1000,50000,1000000)
for (i in 1:length(K)) {
  for(j in 1: length(N_Simul)){
    f.call <- function(x) max(x-K[i], 0)
    f.put <- function(x) max(K[i]-x, 0)
    set.seed(123)
    M <- N_Simul[j]
    MC_Call[i,j] <- MCPrice(x=S0, t=0, T=T, r=r, sigma, M=M, f=f.call)
    MC_Put[i,j] <- MCPrice(x=S0, t=0, T=T, r=r, sigma, M=M, f=f.put)
  }
}

