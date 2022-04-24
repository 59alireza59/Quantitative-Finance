# BLACK-SCHOLES MODELING OF FINANCIAL TIME SERIES

################################################################################
## This task is an attempt to present a framework to check the validity of 
## the assumptions of Black-Scholes model on sample components 
## (ten Euro Area banks) over two different time horizons: 
## a long-term one (ten years of daily stock prices and log-returns) 
## and a short-term one (one year finishing on the 30th April 2015):
##
## Deutsche Bank (DBK.DE); Commerzbank (CBK.DE); BNP Paribas (BNP.PA);
## Credit Agricole (ACA.PA); Banco Santander (SAN.MC); Banco Bilbao (BBVA.MC);
## Intesa San Paolo (ISP.MI); Unicredit (UCG.MI); Mediobanca (MB.MI); and
## ING Groep (INGA.AS).
################################################################################

library(plyr)
library(quantmod)
require(e1071)
library(tseries)
library(TTR)
library(fGarch)


## Check Normality of Log-returns using Geometric Brownian Motion concept:

################################################################################
## Check graphically the distribution of assets 
## by plotting the empirical densities of assets and 
## compare them with the normal distribution 
## with respect to the mean and the standard deviation of log-returns.
################################################################################

## Long-run
## create environment to load data into
data.env <- new.env()

## Access to Market Value Release at Yahoo Finance 
## for multiple assets at one glance:
stocks <- c("ISP.MI", "UCG.MI", "MB.MI", "SAN.MC", "BBVA.MC", 
            "DBK.DE", "CBK.DE", "BNP.PA", "ACA.PA", "INGA.AS")
l_ply(stocks, function(sym) try(getSymbols(sym, from="2005-04-30", to="2015-04-30", 
                                           env=data.env),silent=T))

## Drop all "bad" tickers and merge the good stocks by collecting 
## all the good stock xts objects
stocks <- stocks[stocks %in% ls(data.env)]
market.data <- xts()
for(i in seq_along(stocks)) {
  symbol <- stocks[i]
  market.data <- merge(market.data, Ad(get(symbol,envir=data.env)))
}
colnames(market.data) <- c("ISP.MI", "UCG.MI", "MB.MI", "SAN.MC", "BBVA.MC", 
                           "DBK.DE", "CBK.DE", "BNP.PA", "ACA.PA", "INGA.AS")

## Face to missing data
market.data<-na.omit(diff(log(market.data)))

## Compute the densities of the market data
dens <- apply(market.data, 2, density)

## For each asset display how the normal distribution could or
## could not fit the distribution of the log-returns over the long term
par(mfrow=c(2,5))
for (i in 1:10) {
  plot(dens[[i]], main=names(dens)[i])
  curve(dnorm(x, mean(market.data[,i]), sd(market.data[,i])), add=TRUE, col="red")
}

## check the symmetric status of the distributions and existence of 
## leptokurtic behavior in each case via Kurtosis and Skewness tests

## Skewness test
skewness_data <- apply(market.data, 2, skewness)
print(skewness_data)

## Kurtosis test
kurtosis_data <- apply(market.data, 2, kurtosis)
print(kurtosis_data)

## Normality Assessment by estimating 
## the distribution parameters of the fitting process 
par(mfrow=c(2,5))
for(i in 1:dim(market.data)[2]){
  qqnorm(market.data, xlab=names(market.data)[i])
  qqline(market.data, col="red")
}

## Examine the hypothesis through the p-values of the normality assessment
shapiro <- apply(market.data, 2, shapiro.test)
shapiro


################################################################################
### Short-run
#
## Data Preparation
#stocks <- c("ISP.MI", "UCG.MI", "MB.MI", "SAN.MC", "BBVA.MC", 
#            "DBK.DE", "CBK.DE", "BNP.PA", "ACA.PA", "INGA.AS")
#l_ply(stocks, function(sym) try(getSymbols(sym, from="2005-04-30", to="2015-04-30", 
#                                           env=data.env),silent=T))
#stocks <- stocks[stocks %in% ls(data.env)]
#market.data <- xts()
#for(i in seq_along(stocks)) {
#  symbol <- stocks[i]
#  market.data <- merge(market.data, Ad(get(symbol,envir=data.env)))
#}
#colnames(market.data) <- c("ISP.MI", "UCG.MI", "MB.MI", "SAN.MC", "BBVA.MC", 
#                           "DBK.DE", "CBK.DE", "BNP.PA", "ACA.PA", "INGA.AS")
#market.data<-na.omit(diff(log(market.data)))
#
### Compute densities
#dens <- apply(market.data, 2, density)
#par(mfrow=c(2,5))
#for (i in 1:10) {
#  plot(dens[[i]], main=names(dens)[i])
#  curve(dnorm(x, mean(market.data[,i]), sd(market.data[,i])), add=TRUE, col="red")
#}
#
### Skewness test
#skewness_data <- apply(market.data, 2, skewness)
#print(skewness_data)
#
### Kurtosis test
#kurtosis_data <- apply(market.data, 2, kurtosis)
#print(kurtosis_data)
#
### Normality Assessment
#par(mfrow=c(2,5))
#for(i in 1:dim(market.data)[2]){
#  qqnorm(market.data, xlab=names(market.data)[i])
#  qqline(market.data, col="red")
#}
#
### Sharipo test
#shapiro <- apply(market.data, 2, shapiro.test)
#shapiro
################################################################################

## CHECKING CONSTANT VARIANCE OF LOG-RETURNS per Long-run

################################################################################
## Black and Scholes model assumes variance of the returns for 
## the stocks may constant over time. However, this condition 
## usually doesn’t hold in the reality. 
################################################################################

## Fix the period of long run
start <- "2005-04-30"
end <- "2015-04-30"

## Get the historical data
ISP <- get.hist.quote("ISP.MI",quote="Adjusted",start=start, end=end)
UCG <- get.hist.quote("UCG.MI",quote="Adjusted",start=start, end=end)
MB <- get.hist.quote("MB.MI",quote="Adjusted",start=start, end=end)
SAN <- get.hist.quote("SAN.MC",quote="Adjusted",start=start, end=end)
BBVA <- get.hist.quote("BBVA.MC",quote="Adjusted",start=start, end=end)
DBK <- get.hist.quote("DBK.DE",quote="Adjusted",start=start, end=end)
CBK <- get.hist.quote("CBK.DE",quote="Adjusted",start=start, end=end)
BNP <- get.hist.quote("BNP.PA",quote="Adjusted",start=start, end=end)
ACA <- get.hist.quote("ACA.PA",quote="Adjusted",start=start, end=end)
INGA <- get.hist.quote("INGA.AS",quote="Adjusted",start=start, end=end)
Series <- zoo(cbind(ISP, UCG, MB, SAN, BBVA, DBK, CBK, BNP, ACA, INGA))
colnames(Series) <- gsub("Adjusted.","",colnames(Series))

## Remove NAs
Series<-na.approx(Series,rule=2)
Series_ROC<-ROC(Series)
Series_ROC<-na.approx(Series_ROC,rule=2)

## Plot the log-returns for the fixed period
plot(Series_ROC, xlab="Time", main="LOG-RETURNS FOR 2005-2015 TIME SERIES")

## Test the autocorrelation of log-returns using a Box-Ljung test 
## (e.g. up to the 10th lag) in which Null Hypothesis is the absence of 
## autocorrelation among the elements of a time series.

Box.test_data10 <- apply(Series_ROC, 2, function(x) Box.test(x,lag=10))
Box.test_data10

################################################################################
## Perform Box-Ljung autocorrelation test for other lags as well:
Box.test_data3 <- apply(Series_ROC, 2, function(x) Box.test(x,lag=3))
Box.test_data3
Box.test_data6 <- apply(Series_ROC, 2, function(x) Box.test(x,lag=6))
Box.test_data6
Box.test_data12 <- apply(Series_ROC, 2, function(x) Box.test(x,lag=12))
Box.test_data12
              ###############################################
# Source: https://github.com/cran/FinTS/blob/master/R/autocorTest.R
AutocorTest <- function(x, lag=ceiling(log(length(x))),
                        type=c("Ljung-Box", "Box-Pierce", "rank"),
                        df=lag ){
  type <- match.arg(type)
  if(type=="rank"){
    x <- rank(x) 
    type <- "Ljung-Box"
  }
  #
  LjB <- Box.test(x, lag, type)
  #
  if(df<=0) df <- 1
  if(df != lag){
    LjB$parameter <- c(df=df) 
    #   Correct the p.value  
    LjB$p.value <- pchisq(LjB$statistic, df, lower.tail=FALSE) 
    LjB$method <- paste(LjB$method, " (lag = ", lag, ")", sep="")
  }
  LjB$data.name <- deparse(substitute(x))
  LjB$Total.observ <- length(x)
  #  
  LjB  
}
             ###############################################
My_AutocorTest1 <- apply(Series_ROC, 2, AutocorTest)
My_AutocorTest1

My_AutocorTest2 <- apply(Series_ROC^2, 2, AutocorTest)
My_AutocorTest2

My_AutocorTest3 <- apply(abs(Series_ROC), 2, AutocorTest)
My_AutocorTest3
################################################################################

## Existence of (strong) autocorrelation of returns for time series is 
## an (strong) evidence for volatility clustering and GARCH effect.  

Garch_data11 <- apply(Series_ROC, 2, function(x) garchFit(~garch(1,1), data=x , trace=FALSE))
Garch_data11




