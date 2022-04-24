# MONTE CARLO METHOD AND MODEL SELECTION

## This task use Monte Carlo Method for option pricing
## for American Express Company (AXP) provided by Yahoo Finance.
## Basically, this method generates a random sample for 
## the underlying assets in a risk neutral world by obtaining 
## a set of different payoffs at time in order to calculate 
## the mean of the sample finding an expected payoff.
## The advantage of such utilization is that any type of 
## stochastic processes like GMB or other models can be implemented
## inside the Monte Carlo Model.

## install.packages("pacman")
## Load required packages at once
library(pacman)
pacman::p_load(quantmod, yuima, fOptions,sde)

######################## MCP AND MODEL SELECTION
###### MONTECARLO SIMULATION #####
## GBM SIMULATION##
MCPrice <- function(model, S0, r, T, f, par, M=10000){ set.seed(123)
  samp <- setSampling(Terminal=T, n=100) 
  payoff <- NULL
  for(i in 1:M){
    X <- simulate(model, xinit=S0, sampling=samp, true.par=par)
    y <- get.zoo.data(X)[[1]] 
    # extract the data
    ST <- as.numeric( y[length(y)] )  
    # terminal value
    payoff <- c(payoff, f(ST))
  }
  p <- mean(payoff)*exp(-r*T)
  return(p)
}
f.call <- function(x) max(x-K, 0)
f.put <- function(x) max(K-x, 0)
Delta <- 1/252
S0 <- 105.98   ## the last price of EXP on 2018-09-04 1:17PM EDT
T<- 3/252 ## The option expired on the 07/09
r<- 0.0019 ## 1 year USA t-bill rate
sigma.hat <- 0.2423211 ## historical volatility
getSymbols("AXP", from="2018-02-01", to="2018-09-04")
S <- AXP$AXP.Adjusted
gbm <- setModel(drift="mu*x", diff="sigma*x", solve.variable = c("x"))
yuima <- setYuima(model=gbm, data=setData(S, delta=Delta))
fit.gbm <- qmle(yuima, start=list(mu=1, sigma=0.5))

#################CALL
Ks <- c(102, 103, 104, 105)
#strike prices
nK <-length(Ks)
Ps <- numeric(nK)
#pricing results
for(j in 1 : nK) {
  K <- Ks[j]
  P <- MCPrice(gbm, S0, r, T, f.call, par=as.list(coef(fit.gbm)))
  Ps[j] <- P
}
Ks
Ps
# Ks [1] 102 103 104 105
# Ps [1] 4.233298 3.334054 2.514833 1.800493

################ PUT
Ks <- c(102, 103, 104, 105)
#strike prices
nK <-length(Ks)
Ps <- numeric(nK)
#pricing results
for(j in 1 : nK) {
  K <- Ks[j]
  P <- MCPrice(gbm, S0, r, T, f.put, par=as.list(coef(fit.gbm)))
  Ps[j] <- P
}
Ks
Ps
# Ks [1] 102 103 104 105
# Ps [1] 0.08901315 0.18974656 0.37050241 0.65614034


################# MODEL SELECTION (AFTER CHANGING POINT)
getSymbols("AXP", from="2018-02-01", to="2018-09-04")
S <- AXP$AXP.Adjusted
S
S <- AXP$AXP.Close
S
X <- na.omit(diff(log(S)))
m <- mean(S) 
m
#  98.61791
plot(S,main="AXP STOCK PRICE")
abline(h = m, lty = 1, col = "red", lwd=1.5) 
axis(4, m, expression(mean))
plot(X) # tesla log returns
s <- as.numeric(sqrt(var(X)*252)) 
s
#  0.2423211
M <- mean(X)*252 + 0.5*s^2 
M
#  0.1289259
## GBM ##
gbm <- setModel(drift="mu*x", diff="sigma*x", solve.variable = c("x"))
yuima <- setYuima(model=gbm, data=setData(S,delta=Delta))
upper <- list (mu=0.5, sigma=0.6)
lower <-list (mu=0.2, sigma=0.01)
# gbm.fit <- qmle(yuima, start=list(mu=1, sigma=1))
fit.gbm <- qmle(yuima, start=list(mu=M, sigma=s),upper=upper,lower=lower,method="L-BFGS-B")
summary(fit.gbm)
## VASICEK ##
vas <- setModel(drift="alpha*(mu-x)", diff="sigma", solve.variable = c("x"))
yuima.vas <- setYuima(model=vas, data=setData(S, delta=Delta))
#alpha= speed of convergence, compresa tra 0 e 1
#upper <- list (mu=100, sigma=0.6, alpha=0.5)
#lower <-list (mu=90, sigma=0.01, alpha=0.5)
fit.vas <- qmle(yuima.vas, start=list(alpha=0.5, mu=m, sigma=s),method="L-BFGS-B")
summary(fit.vas)
## CIR ##
cir <- setModel(drift="alpha*(mu-x)", diff="sigma*sqrt(x)", solve.variable = c("x"))
yuima.cir <- setYuima(model=cir, data=setData(S, delta=Delta))
fit.cir <- qmle(yuima.cir, start=list(alpha=0.5, mu=m, sigma=s),method="L-BFGS-B")
summary(fit.cir)
## CKLS ##
ckls <- setModel(drift="alpha*(mu-x)", diff="sigma*x^beta", solve.variable = c("x"))
yuima.ckls <- setYuima(model=ckls, data=setData(S, delta=Delta))
fit.ckls <- qmle(yuima.ckls, start=list(alpha=0.5, mu=m, sigma=s, beta=0.5),method="L-BFGS-B")
summary(fit.ckls)
## GBM1 ##  WITH MEAN REVERSION 
#variaton of GBM with the drift component from VAS
gbm1 <- setModel(drift="mu-x", diff="sigma*x", solve.variable = c("x"))
yuima.gbm1 <- setYuima(model=gbm1, data=setData(S, delta=1/252))
fit.gbm1 <- qmle(yuima.gbm1, start=list(mu=m, sigma=s))
summary(fit.gbm1)
AIC(fit.gbm)  #  540.2773
AIC(fit.vas)  #  532.2503
AIC(fit.cir)  #  534.8684
AIC(fit.ckls) #  529.4595
AIC(fit.gbm1) #  539.6496
#according to the AIC, the CKLS is the best model.

################ CKLS SIMULATION##
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
f.call <- function(x) max(x-K, 0)
f.put <- function(x) max(K-x, 0)
Delta <- 1/252
S0 <- 105.98   ## the last price of EXP on 2018-09-04 1:17PM EDT
T<- 3/252 ## The option expired on the 07/09
r<- 0.0019 ## 1 year USA t-bill rate
sigma.hat <- 0.2423211 # historical volatility
getSymbols("AXP", from="2018-02-01", to="2018-09-04")
S <- AXP$AXP.Adjusted
ckls <- setModel(drift="alpha*(mu-x)", diff="sigma*x^beta", solve.variable = c("x"))
yuima.ckls <- setYuima(model=ckls, data=setData(S, delta=Delta))
fit.ckls <- qmle(yuima.ckls, start=list(alpha=0.5, mu=m, sigma=s, beta=0.5),method="L-BFGS-B")
#################CALL
Ks <- c(102, 103, 104, 105)
#strike prices
nK <-length(Ks)
Pckls.c <- numeric(nK)
#pricing results
for(j in 1 : nK) {
  K <- Ks[j]
  P <- MCPrice(ckls, S0, r, T, f.call, par=as.list(coef(fit.ckls)))
  Pckls.c[j] <- P
}
Ks
Pckls.c
# Ks [1]  102 103 104 105
# Pckls.c [1] 3.0833214 2.2035632 1.4411584 0.8409291
################PUT
Ks <- c(102, 103, 104, 105)
#strike prices
nK <-length(Ks)
Pckls.p <- numeric(nK)
#pricing results
for(j in 1 : nK) {
  K <- Ks[j]
  P <- MCPrice(ckls, S0, r, T, f.put, par=as.list(coef(fit.ckls)))
  Pckls.p[j] <- P
}
Ks
Pckls.p
# Ks [1]  102 103 104 105
# Pckls.p [1]  0.07818804 0.19840715 0.43597978 0.83572784
