# MONTE CARLO METHOD AND MODEL SELECTION

################################################################################
## This task use Monte Carlo Method for option pricing
## for American Express Company (AXP) provided by Yahoo Finance.
## Basically, this method generates a random sample for 
## the underlying assets in a risk neutral world by obtaining 
## a set of different payoffs at time in order to calculate 
## the mean of the sample finding an expected payoff.
## The advantage of such utilization is that any type of 
## stochastic processes like GMB or other models can be implemented
## inside the Monte Carlo Model.
################################################################################

## install.packages("pacman")
## Load required packages at once
library(pacman)
pacman::p_load(quantmod, yuima, fOptions,sde)


## Simulate Geometric Brownian motion using Monte Carlo method and model selection

## Build Monte Carlo Price
MCPrice <- function(model, S0, r, T, f, par, M=10000){ set.seed(123)
  samp <- setSampling(Terminal=T, n=100) 
  payoff <- NULL
  for(i in 1:M){
    X <- simulate(model, xinit=S0, sampling=samp, true.par=par)
    y <- get.zoo.data(X)[[1]] 
    # extract data
    ST <- as.numeric( y[length(y)] )  
    # terminal value
    payoff <- c(payoff, f(ST))
  }
  p <- mean(payoff)*exp(-r*T)
  return(p)
}

## Generate Parameters of Geometric Brownian motion
f.call <- function(x) max(x-K, 0)
f.put <- function(x) max(K-x, 0)
Delta <- 1/252
S0 <- 105.98   ## Last price of EXP on 2018-09-04 1:17PM EDT
T<- 3/252 ## Option expired on the 07/09
r<- 0.0019 ## 1 year USA t-bill rate
sigma.hat <- 0.2423211 ## historical volatility

## Load stock data from yahoo finance
getSymbols("AXP", from="2018-02-01", to="2018-09-04")
S <- AXP$AXP.Adjusted

## Simulate Geometric Brownian motion
gbm <- setModel(drift="mu*x", diff="sigma*x", solve.variable = c("x"))
yuima <- setYuima(model=gbm, data=setData(S, delta=Delta))
fit.gbm <- qmle(yuima, start=list(mu=1, sigma=0.5))

## Get Call option
Ks <- c(102, 103, 104, 105)

## Get strike prices
nK <-length(Ks)
Ps <- numeric(nK)

## Obtain pricing results
for(j in 1 : nK) {
  K <- Ks[j]
  P <- MCPrice(gbm, S0, r, T, f.call, par=as.list(coef(fit.gbm)))
  Ps[j] <- P
}
Ks
Ps
# Ks [1] 102 103 104 105
# Ps [1] 4.233298 3.334054 2.514833 1.800493

## Get Put option 
Ks <- c(102, 103, 104, 105)

## Get strike prices
nK <-length(Ks)
Ps <- numeric(nK)

## Obtain pricing results
for(j in 1 : nK) {
  K <- Ks[j]
  P <- MCPrice(gbm, S0, r, T, f.put, par=as.list(coef(fit.gbm)))
  Ps[j] <- P
}
Ks
Ps
# Ks [1] 102 103 104 105
# Ps [1] 0.08901315 0.18974656 0.37050241 0.65614034


## Perform model selection procedure
## Load stock data
getSymbols("AXP", from="2018-02-01", to="2018-09-04")
S <- AXP$AXP.Adjusted
S
S <- AXP$AXP.Close
S
X <- na.omit(diff(log(S)))
m <- mean(S) 
m #  98.61791

## Plot historical data
plot(S,main="AXP STOCK PRICE")
abline(h = m, lty = 1, col = "red", lwd=1.5) 
axis(4, m, expression(mean))
plot(X) # AXP log returns

## See brief descriptive statistic
s <- as.numeric(sqrt(var(X)*252)) 
s #  0.2423211
M <- mean(X)*252 + 0.5*s^2 
M #  0.1289259


## GBM construction
gbm <- setModel(drift="mu*x", diff="sigma*x", solve.variable = c("x"))
yuima <- setYuima(model=gbm, data=setData(S,delta=Delta))
upper <- list (mu=0.5, sigma=0.6)
lower <-list (mu=0.2, sigma=0.01)
# gbm.fit <- qmle(yuima, start=list(mu=1, sigma=1))
fit.gbm <- qmle(yuima, start=list(mu=M, sigma=s),upper=upper,lower=lower,method="L-BFGS-B")
summary(fit.gbm)

## Create VASICEK method
vas <- setModel(drift="alpha*(mu-x)", diff="sigma", solve.variable = c("x"))
yuima.vas <- setYuima(model=vas, data=setData(S, delta=Delta))
#alpha= speed of convergence, compresa tra 0 e 1
#upper <- list (mu=100, sigma=0.6, alpha=0.5)
#lower <-list (mu=90, sigma=0.01, alpha=0.5)
fit.vas <- qmle(yuima.vas, start=list(alpha=0.5, mu=m, sigma=s),method="L-BFGS-B")
summary(fit.vas)

## Create CIR method 
cir <- setModel(drift="alpha*(mu-x)", diff="sigma*sqrt(x)", solve.variable = c("x"))
yuima.cir <- setYuima(model=cir, data=setData(S, delta=Delta))
fit.cir <- qmle(yuima.cir, start=list(alpha=0.5, mu=m, sigma=s),method="L-BFGS-B")
summary(fit.cir)

## Create CKLS method
ckls <- setModel(drift="alpha*(mu-x)", diff="sigma*x^beta", solve.variable = c("x"))
yuima.ckls <- setYuima(model=ckls, data=setData(S, delta=Delta))
fit.ckls <- qmle(yuima.ckls, start=list(alpha=0.5, mu=m, sigma=s, beta=0.5),method="L-BFGS-B")
summary(fit.ckls)

## Perform GBM simulation with mean reversion
gbm1 <- setModel(drift="mu-x", diff="sigma*x", solve.variable = c("x"))
yuima.gbm1 <- setYuima(model=gbm1, data=setData(S, delta=1/252))
fit.gbm1 <- qmle(yuima.gbm1, start=list(mu=m, sigma=s))
summary(fit.gbm1)
AIC(fit.gbm)  #  540.2773
AIC(fit.vas)  #  532.2503
AIC(fit.cir)  #  534.8684
AIC(fit.ckls) #  529.4595, CKLS is the best model 
AIC(fit.gbm1) #  539.6496

## CKLS simulation
## Build funaction of Make Monte Carlo pricing
MCPrice <- function(model, S0, r, T, f, par, M=10000){ set.seed(123)
  samp <- setSampling(Terminal=T, n=100) 
  payoff <- NULL
  for(i in 1:M){
    X <- simulate(model, xinit=S0, sampling=samp, true.par=par)
    y <- get.zoo.data(X)[[1]] 
    ST <-  as.numeric( y[length(y)] ) 
    payoff <- c(payoff, f(ST))
  }
  p <- mean(payoff)*exp(-r*T)
  return(p)
}

## Get call and put options
f.call <- function(x) max(x-K, 0)
f.put <- function(x) max(K-x, 0)
Delta <- 1/252
S0 <- 105.98   ## Last price of EXP on 2018-09-04 1:17PM EDT
T<- 3/252 ## Option expired on the 07/09
r<- 0.0019 ## 1 year USA t-bill rate
sigma.hat <- 0.2423211 # historical volatility
getSymbols("AXP", from="2018-02-01", to="2018-09-04")
S <- AXP$AXP.Adjusted
ckls <- setModel(drift="alpha*(mu-x)", diff="sigma*x^beta", solve.variable = c("x"))
yuima.ckls <- setYuima(model=ckls, data=setData(S, delta=Delta))
fit.ckls <- qmle(yuima.ckls, start=list(alpha=0.5, mu=m, sigma=s, beta=0.5),method="L-BFGS-B")

## Generate Call options
Ks <- c(102, 103, 104, 105)

## Simulate strike prices
nK <-length(Ks)
Pckls.c <- numeric(nK)

## Get pricing results
for(j in 1 : nK) {
  K <- Ks[j]
  P <- MCPrice(ckls, S0, r, T, f.call, par=as.list(coef(fit.ckls)))
  Pckls.c[j] <- P
}
Ks
Pckls.c
# Ks [1]  102 103 104 105
# Pckls.c [1] 3.0833214 2.2035632 1.4411584 0.8409291

## Generate put option
Ks <- c(102, 103, 104, 105)

## Simulate strike prices
nK <-length(Ks)
Pckls.p <- numeric(nK)

## Get pricing results
for(j in 1 : nK) {
  K <- Ks[j]
  P <- MCPrice(ckls, S0, r, T, f.put, par=as.list(coef(fit.ckls)))
  Pckls.p[j] <- P
}
Ks
Pckls.p
# Ks [1]  102 103 104 105
# Pckls.p [1]  0.07818804 0.19840715 0.43597978 0.83572784
