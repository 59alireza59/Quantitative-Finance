# LEVY MARKET PROCESS

################################################################################
## Unlike the unrealistic assumptions used in the Black and Scholes model - 
## the particular stochastic process (GBM) which lies below the stock prices, 
## the normal distribution of the log-returns and the constant volatility - 
## we consider a different method so-called Lévy process where our distributions
## are not normal. Generally, Lévy process is defined by a family of processes
## wide enough to comprise a variety of well-known other stochastic processes.
## As a stochastic process, it can be defined by the sum of a jump process
## (the Poisson process) and a Brownian motion with drift characterizing by
## similar properties to those of the Brownian motion, with the relevant
## difference in the distribution of the increments, which are no more Gaussian. 
################################################################################

## install.packages("pacman")
## Load required packages at once
library(pacman)
pacman::p_load(quantmod, yuima, tseries,fBasics)


## Simulate Levy Process

## Access stock data from yahoo finance
getSymbols("AXP", from="2018-02-01", to="2018-09-04")
S <- AXP$AXP.Close
X <- returns(S)

## See AXP prices and log-returns paths
par(mfrow=c(1,2))
lineChart(S,layout=NULL,theme="white",name="AXP Prices")
lineChart(X,layout=NULL,theme="white",name="AXP Log-Returns")
X <- X[!is.na(X)]
sigma.hat <- sqrt( var(X)/(deltat(S)) )
alpha.hat <- mean(X)/(deltat(S))
mu.hat <- alpha.hat + 0.5 * sigma.hat^2
sigma.hat # 0.01526479
mu.hat # 0.0005116106

## Plot the AXP returns density distribution versus the Gaussian distribution
plot(density(X))
d <- function(x) dnorm(x, mean=mu.hat, sd=sigma.hat)
curve(d, -.2, .2, col="red",add=TRUE, n=500)

## Plot the AXP returns density distribution versus different distributions
par(mfrow=c(2,2))
par(mar=c(3,3,3,1))
grid <- NULL

# Fit different distributions to stock data
nFit(X)
nigFit(X,trace=FALSE)
hypFit(X,trace=FALSE)
tFit(X,trace=FALSE)
#ghFit(X,trace=FALSE)
nigFit(X,trace=FALSE)

## Plot the AXP returns density distribution versus the NIG distribution:
## Load data and Prepare data
getSymbols("AXP", from="2018-02-01", to="2018-09-04")
S <- AXP$AXP.Close
X <- returns(S)
S <- rev(as.numeric(S[,1]))
X <- returns(S)
n <- length(S)
n # 148
S1 <- S[1:(n*0.85)]
# Split time series
m <- length(S1)
m # 125
S.true <- S[-(1:(m-1))]
S.true
N <- length(S.true)
N # 24
X <- as.numeric(returns(S1))
X <- X[-1]
# Fit the model 
X <- as.numeric(X)
X
mod.NIG <- nigFit(X,trace=FALSE)@fit$par # Consider the fact that mod.NIG is a vector of number
alpha   <- mod.NIG["alpha"]
beta    <- mod.NIG["beta"]
delta   <- mod.NIG["delta"]
mu      <- mod.NIG["mu"]
xvals <- seq(from=min(X),to=max(X),length=100);
yvals.NIG <- dnig(xvals,alpha=alpha,beta=beta,delta=delta,mu=mu)
plot(density(X))  # black= distr of returns
lines(xvals,yvals.NIG,col="red",lty=5,lwd=2) # red= distr of NIG
lines(xvals, dnorm(xvals, mean=mean(X), sd=sd(X)), col="blue", lty=3,lwd=2) # blue= GBM distribution
S0 <- S1[m] # It is terminal value of the first part of data
S0 # 95.64

## Build Levy process code
obj.model <- setModel(drift="0", diff="0", jump.coeff="1", measure.type="code", 
                      measure=list(df=sprintf("rNIG(z, alpha=%f, beta=%f, delta=%f, mu=%f, Lambda=Lambda )",
                                              alpha,beta,delta,mu)))
Lambda <- matrix(1,1,1)
obj.sampling <- setSampling(Terminal=N, n=N)
obj.yuima <- setYuima(model=obj.model, sampling=obj.sampling)
result <- simulate(obj.yuima,xinit=0, true.par = list(Lambda=Lambda))
Xt <- result@data@zoo.data[[1]]
SS <- S0*exp(Xt)
## Determine the trajectory of the model
A <- c(S1, S.true)
plot(A,type="n") # AXP simulated pattern under Lévy process (red line) and real pattern (blue line).
points(S1,type="l")
points(m+1:length(SS) , SS,type="l",col="red") # simulated data (trajectory)
points(m+1:length(S.true), S.true,type="l",col="blue") #true data of the model

## Simulate process of pricing with Levy
set.seed(123)
M <- 1000
sims <- numeric(M)
for(i in 1:M){
  result <- simulate(obj.yuima,xinit=0,true.par=list(Lambda=Lambda))
  Xt <- result@data@zoo.data[[1]]
  SS <- S0*exp(Xt)
  sims[i] <- SS[length(SS)]
}
A <- c(S1, S.true)
A
plot(A,type="n",ylim=c(180,280))
points(S1,type="l")
points(m+1:length(SS) , SS,type="l",col="red")
points(m+1:length(S.true), S.true,type="l",col="blue")
points(rep(length(A),M), sims) #all the points are simulated values of the model
f<- function(x) max(x-K, 0)

## Make call payoff function
Delta <- 1/252
K <- 105.98
r <- 0.0019
# T <- N-1
T <- 3*Delta
p <- mean(sapply(sims, f))
p # 0.06857888 # payoff under Levy NIG model
p <- mean(sapply(sims, f))*exp(-r*T) #discounting
p # 0.06857733
s <- sqrt(var(X)/Delta)
s # 0.2071822
mu <- mean(X)/Delta + 0.5 * s^2
mu # -0.1871679
MCP <- function(x=1, t=0, T=1, r=1, sigma=1, M=1000, f){ set.seed(123)
  u <- rnorm(M)
  h <-  x*exp((r-0.5*sigma^2)*(T-t)+sigma*sqrt(T-t)*u)
  p <- sapply(h, function(x) f(x))
  p <- mean(p)
  #
  p*exp(-r*(T-t)) 
  p
}
p0 <- MCP(S0, 0, T=T, r=r, sigma=s, M=M, f=f)
p0 # 0 # B&S expected payoff
## Compare these quantities:
levy.price <- sapply(sims, f)

## Plot Lévy payoff
plot(density(levy.price),xlim=c(-5,10), main="Density Levy Price") # Lévy payoff (red line), GBM payoff (blue line) with respect to the null real payoff.
abline(v=f(S[length(S)])) #black true payoff trajectory
abline(v=p,col="red") #red levy payoff
abline(v=p0,col="blue") #blue B&S payoff
f(S[length(S)])  #0 real payoff

