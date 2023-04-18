# Quantitative-Finance
This Quantitative Finance Engineering project involves option pricing and estimation of financial models using R

## Abstract:
This project programms an application of option pricing theory to real data, calculating various operations on stock and option prices. Despite the existence of some unacceptable assumptions in reality, which lead to some weaknesses in the theoretical interpretation, this programming work attempts instead to provide a valuable statistical framework with an analytical perspective. To this end, various exploratory techniques, such as clustering and lead-laggard estimation, have been used to identify some similarities between well-known firms in the financial services market and, accordingly, to show which of them is likely to be the market leader. The program also checks whether and how the values obtained by using the Black and Scholes model and the Monte Carlo method to value some options differ from the corresponding market prices. These values also allow us to identify the main reasons for this divergence and the assumptions on which they are based. The AIC method is then used to find the best model for interpreting the path of the time series of prices, taking into account some well-known models: Geometric Brownian Motion (GBM), Vasicek (VAS), Cox-Ingersoll-Ross (CIR) and CKLS. Finally, another option pricing method, the Lèvy market model, is used. The main reason for choosing the Lèvy market model is that it takes into account all the features that lie outside the standard Black and Scholes assumptions, such as some typical price features such as jumps, changes in volatility and non-normality of returns.

## Intro:

This platform is programmed to demonstrate in practice which statistical concept can be developed for which part of financial market analysis:

1- Time Series Clustering

This task propounds a clustering framework to draw the main features of the Johnson & Johnson’s stock exchange compared with its main competitors. To build the framework, we considered a time series of daily closing prices from Yahoo Finance for two different sectors (pharmaceutical/healthcare/biotechnology sector and banking services, from September 2017 to September 2018 for the following financial assets):
- Pfizer (PFE);
- Novartis (NVS);
- Teva Pharmaceutical Industries Limited (TEVA);
- Amgen Inc. (AMGN);
- Sanofi (SNY);
- GlaxoSmithKline (GSK);
- Merck & Co. (MRK);
- Viatris Inc. (VTRS);
- JP Morgan (JPM);
- HSBC Holding (HSBC); and
- Nomura (NMR).

2- Lead–lag Effect

The lead-lag effect illustrates how the lagging price process tends to mimic the fluctuations of the leading price process with a temporary lag, or vice versa. The aim of this assignment is to identify the main characteristics between an American global investment management company called BlackRock and its competitors by calculating correlation matrices and lead-lag estimation.

3- Black–Scholes Test

This exercise is an attempt to provide a framework for checking the validity of the assumptions of the Black-Scholes model on sample components (ten euro area banks) over two different time horizons - a long-term one (ten years of daily share prices and logarithmic returns) and a short-term one (one year ending 30 April 2015):
- Deutsche Bank (DBK.DE);
- Commerzbank (CBK.DE);
- BNP Paribas (BNP.PA);
- Credit Agricole (ACA.PA);
- Banco Santander (SAN.MC);
- Banco Bilbao (BBVA.MC);
- Intesa San Paolo (ISP.MI);
- Unicredit (UCG.MI);
- Mediobanca (MB.MI); and
- ING Groep (INGA.AS).

4- Change Point Detection

This exercise attempts to determine "whether or not the stock pattern of Wells Fargo follows a geometric Brownian motion" by detecting changes. In answering this question, both the independence and the normality of the logarithmic returns must be checked in the presence of the Black & Scholes assumption.

5- Option Pricing

This exercise presents a framework for option pricing using the Black and Scholes formula and Monte Carlo simulation, where the first requirement is to apply Changing Point Analysis to a rich financial time series dataset from Yahoo Finance (Johnson & Johnson stock prices). The Changing Point Analysis helps to identify when there is the volatility changes in the data relative to the length of the desired duration.

6- Monte Carlo Method

This exercise uses the Monte Carlo method to calculate option prices for American Express Company (AXP), provided by Yahoo Finance. Basically, this method generates a random sample for the underlying assets in a risk-neutral world by obtaining a set of different payoffs at each point in time, and then calculating the average of the sample that finds an expected payoff. The advantage of this is that any type of stochastic process, such as GMB or other models, can be implemented within the Monte Carlo model.

7- Levy Market Process

In contrast to the unrealistic assumptions used in the Black & Scholes model - the particular stochastic process (GBM) underlying stock prices, the normal distribution of the log-returns and the constant volatility - this extract considers a an alternative method, the so-called Lévy process where the distributions are not normal. In general, the Lévy process is defined by a family of processes broad enough to include a variety of other well-known stochastic processes. As a stochastic process, it can be defined as the sum of a jump process (the Poisson process) and a Brownian motion with drift, characterised by similar properties to those of the Brownian motion, with the relevant difference that the distribution of the increments is no longer Gaussian. 
