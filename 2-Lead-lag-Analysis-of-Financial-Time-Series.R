# LEAD-LAG ANALYSIS OF FINANCIAL TIME SERIES

## install.packages("pacman")
## Load required packages at once
library(pacman)
pacman::p_load(quantmod,plyr,yuima,corrplot)

## The lead–lag effect illustrates how the lagger price process 
## tends to emulate the oscillations of the leader price process 
## relatively with some temporary delay, or vice versa. 
## This task aims at exposing the main features between 
## an American global investment management corporation so-called 
## BlackRock and its competitors through correlation matrix calculation 
## and the lead-lag estimation.

## create environment to load data
data.env <- new.env()


## Load data from Yahoo Finance and Show the historical stock chart
getSymbols("BLK", from="2014-01-03", to="2018-08-31")
chart_Series(Cl(BLK))


## Access to Market Value Release for multiple assets at one glance:
stocks <- c("BLK", "KKR", "XOM", "BP", "F", "CVX", "FB", 
            "TTM","GM", "T","VZ","PEUG.PA", "TM", "STT")
l_ply(stocks, function(sym) try(getSymbols(sym,env=data.env),silent=T))

## Drop all "bad" tickers
stocks <- stocks[stocks %in% ls(data.env)]

## Merge our good stocks by collecting all the good stock xts() objects
market.data <- xts()
for(i in seq_along(stocks)) {
  symbol <- stocks[i]
  market.data <- merge(market.data, Ad(get(symbol,envir=data.env)))
}
colnames(market.data) <- c("BLK", "KKR", "XOM", "BP", "F", "CVX", "FB", 
                           "TTM","GM", "T","VZ","PEUG.PA", "TM", "STT")

## Release statistical relationships among market's participants
## by creating their correlation matrix
mkt<-setYuima(data=setData(market.data,delta=1/252))
round(cce(mkt)$cormat,2)
cols <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow",
                           "white","cyan", "#007FFF","blue","#00007F"))

## Estimate the covariance among the participants
corrplot(cce(mkt)
         $cormat,method="ellipse",
         cl.pos="b",tl.pos="d",tl.srt=60,
         col=cols(100),outline=TRUE)

## Estimate the lead-lag effect in the market
corrplot(llag(mkt),method="ellipse",is.corr=FALSE,
         cl.pos = "b", tl.pos = "d", tl.srt = 60,
         col=cols(100), outline=TRUE)
