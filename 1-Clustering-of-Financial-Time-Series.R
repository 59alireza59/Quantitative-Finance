# CLUSTERING OF FINANCIAL TIME SERIES

## install.packages("pacman")
## Load required packages at once
library(pacman)
pacman::p_load(quantmod,tseries,sde,dtw,proxy,rgl)

## Cluster Analysis of Johnson and Johnson stock:
## This task propounds a clustering framework to 
## draw the main features of the Johnson & Johnson’s 
## stock exchange compared with its main competitors.

## create environment to load data
DataEnv <- new.env()

## To build the framework, a time series of daily closing quotes,
## from yahoo finance for two different sectors: 
## Pharmaceutical/Healthcare/Biotechnology sector and Banking Services,
## have been considered from September 2017 to September 2018
## for the following financial assets:
## Pfizer (PFE); Novartis (NVS);
## Teva Pharmaceutical Industries Limited (TEVA);
## Amgen Inc. (AMGN); Sanofi (SNY);
## GlaxoSmithKline (GSK); Merck & Co. (MRK);
## Viatris Inc. (VTRS); JP Morgan (JPM);
## HSBC Holding (HSBC); and Nomura (NMR).

## Load data and Show the historical stock chart
Symbols <- c("PFE","NVS","TEVA","AMGN","SNY","GSK","MRK",
             "VTRS","JPM","HSBC","NMR")
getSymbols(Symbols, from="2017-09-01", to="2018-09-04", env=DataEnv)    
getSymbols("JNJ", from="2017-09-01", to="2018-09-04")    
chart_Series(Cl(JNJ))
eapply(DataEnv, function(x) add_TA(Cl(x), on = 1))

## Access to Market Value Release:
start <- "2017-09-01"
end <- "2018-09-04"
Johnson <- get.hist.quote("JNJ",quote="Close",start=start, end=end)
Pfizer <- get.hist.quote("PFE",quote="Close",start=start, end=end)
Novartis <- get.hist.quote("NVS",quote="Close",start=start, end=end)
Teva <- get.hist.quote("TEVA",quote="Close",start=start, end=end)
Amgen <- get.hist.quote("AMGN",quote="Close",start=start, end=end)
Sanofi <- get.hist.quote("SNY",quote="Close",start=start, end=end)
GSKline <- get.hist.quote("GSK",quote="Close",start=start, end=end)
Merck <- get.hist.quote("MRK",quote="Close",start=start, end=end)
Viatris <- get.hist.quote("VTRS",quote="Close",start=start, end=end)
JPMorgan <- get.hist.quote("JPM",quote="Close",start=start, end=end)
HSBC <- get.hist.quote("HSBC",quote="Close",start=start, end=end)
Nomura <- get.hist.quote("NMR",quote="Close",start=start, end=end)

Series <- zoo(cbind(Johnson, Pfizer, Novartis, Teva, Amgen, Sanofi, GSKline, Merck, Viatris, JPMorgan, HSBC, Nomura))
colnames(Series) <- c("Johnson", "Pfizer", "Novartis", "Teva", "Amgen", "Sanofi", "GSKline", "Merck", "Viatris", "JPMorgan", "HSBC", "Nomura")

Series <- na.approx(Series, rule=2) ## Face to missing data 

## Plot the price of financial assets to catch the similar trends:
plot(Series, main="Full Set Prices Data", xlab="TIME")

## See the hierarchical structure of Market Players:
## Compute the distance metrics for the hierarchical Clustering
### Markov Operator Distance 
D_MO <- MOdist(Series)
D_MO <- D_MO/max(D_MO)
D_MO
N<-dim(Series)[1]
nSeries<-dim(Series)[2]

### Euclidean Distance
D_initial <- matrix(0,nSeries,nSeries)
for(i in 1:(nSeries - 1))
  for(j in (i+1):nSeries){
    D_initial[i,j] <- sqrt(sum((Series[,i] - Series[,j])^2))
    D_initial[j,i] <- D_initial[i,j]
  }
D_initial <- D_initial/max(D_initial)
colnames(D_initial) <- colnames(Series)
rownames(D_initial) <- colnames(Series)
D_initial
D_E <- matrix(0, nSeries, nSeries)
DELTA <- deltat(Series)

### Short Time Series (STS) Distance
for(i in 1:(nSeries - 1)){
  for(j in (i+1):nSeries){
    D_E[i,j] <- sqrt(sum((diff(Series[,i])/DELTA - diff(Series[,j])/DELTA)^2))
    D_E[j,i] <- D_E[i,j]
  }
}
max_DE <- which.max(D_E)
D_E <- D_E/max_DE
colnames(D_E) <- colnames(Series)
rownames(D_E) <- colnames(Series)
D_E

### Dynamic Time Warping (DTW) Distance
D_DTW <- dist(t(Series),method="dtw")
D_DTW <- D_DTW/max(D_DTW)


### Perform cluster Analysis through dendograms
par(mfrow=c(2,2)); par(mar=c(1,3,3,0))
cl <- hclust(as.dist(D_MO))
plot(cl,main="Markov Operator Distance",xlab="")
cl1 <- hclust(as.dist(D_initial))
plot(cl1, main="Euclidean Distance",xlab="",ylim=c(0,1))
cl2 <- hclust(as.dist(D_E))
plot(cl2, main="STS Distance",xlab="",ylim=c(0,1))
cl3 <- hclust(as.dist(D_DTW))
plot(cl3, main="DTW Distance",xlab="",ylim=c(0,1))

## Plot dendogram of the assets clustered under the Markov Operator distance
clu <- hclust(D_MO)
plot(clu, main="Markov Operator Distance",ylim=c(-1,1))
rect.hclust(cl, k=3, which=c(), border="blue")
gr <- cutree(clu,k=7)
gr

## Do extra job by cluster Analysis in 2/3 Dimensions:

## 2D cluster graph
tmp <- cmdscale(D_MO)
plot(tmp, col=gr, pch="o", main = "Cluster 2D")
text(tmp, rownames(tmp),cex =0.7)

## 3D cluster graph
tmp3<-cmdscale(D_MO,3)
plot3d(tmp3, type="s",col=gr, size=1, cex=2, main="Cluster 3D")
text3d(tmp3,texts=rownames(tmp3), font=5)
