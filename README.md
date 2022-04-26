# Quantitative-Finance
This Quantitative Finance Engineering Project contains Option Pricing and Estimation of Financial Models with R

## Abstract:
This work programs an application of option pricing theory to real data computing different operations on stock and option prices. Despite the existence of some unacceptable assumptions in real leading to some weaknesses for the theory interpretation, this programming work attempts to present a valuable statistical framework with analytical perspective, instead. For doing this, different explorative techniques i.e. clustering and lead-laggards estimation were used to detect some similarities among famous companies in the financial services market, and correspondingly, it shows which of them likely is the market leader as well. Furthermore, this programming check if and how values corresponded to applications of the Black and Scholes model and the Monte Carlo method to price some options, diverge from corresponding market prices. Those values also let us to determine what the main reasons for this divergence are and based on which assumptions. Then, to figure out the best model that to interpret the path followed by the time series of the prices, the AIC method is hired considering some well-known models: Geometric Brownian motion (GBM), Vasicek (VAS), Cox–Ingersoll–Ross (CIR) and CKLS. Finally, a different method for option pricing, so-called the Lèvy market model is used. The main reason to select Lèvy market model is that the model takes into account all characteristics that lie outside of the standard Black and Scholes’ assumptions such as some typical feature of prices as jumps, changes in volatility and non-normality of returns.

## Intro:

This platform is programmed to demonstrate in practice which statistical concept can be developed for which part of financial market analysis:

1- Time Series Clustering

This task propounds a clustering framework to draw the main features of the Johnson & Johnson’s stock exchange compared with its main competitors. To build the framework, a time series of daily closing quotes, from yahoo finance for two different sectors: Pharmaceutical/Healthcare/Biotechnology sector and Banking Services, have been considered from September 2017 to September 2018 for the following financial assets:
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

The lead–lag effect illustrates how the lagger price process tends to emulate the oscillations of the leader price process  relatively with some temporary delay, or vice versa. This task aims at exposing the main features between an American global investment management corporation so-called BlackRock and its competitors through correlation matrix calculation and the lead-lag estimation.

3- Black–Scholes Test

This task is an attempt to present a framework to check the validity of the assumptions of Black-Scholes model on sample components (ten Euro Area banks) over two different time horizons- a long-term one (ten years of daily stock prices and log-returns) and a short-term one (one year finishing on the 30th April 2015):
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

This task is an attempt to check "whether the pattern of the Wells Fargo stock follows a Geometric Brownian Motion or not" due to detecting changes. Answering this question needs to check the both independency and normality of the log returns in the presence of Black & Scholes assumption.

5- Option Pricing

This task represents a framework to compute option price using Black and Scholes formula and Monte Carlo simulation in which the first step requirement is to apply Changing Point Analysis over a rich financial time series dataset from Yahoo Finance (Johnson & Johnson stock prices). The changing point analysis helps to find when there is the volatility changes in data with respect to the length of the desired duration.

6- Monte Carlo Method

This task utilizes the Monte Carlo Method to compute option pricing for American Express Company (AXP) provided by Yahoo Finance. Basically, this method generates a random sample for the underlying assets in a risk neutral world by obtaining a set of different payoffs at time in order to calculate the mean of the sample finding an expected payoff. The advantage of such utilization is that any type of stochastic processes like GMB or other models can be implemented inside the Monte Carlo Model.

7- Levy Market Process

Unlike the unrealistic assumptions used in the Black and Scholes model - the particular stochastic process (GBM) which lies below the stock prices, the normal distribution of the log-returns and the constant volatility - this task considers a different method so-called Lévy process where the distributions are not normal. Generally, Lévy process is defined by a family of processes wide enough to comprise a variety of well-known other stochastic processes. As a stochastic process, it can be defined by the sum of a jump process (the Poisson process) and a Brownian motion with drift characterizing by similar properties to those of the Brownian motion, with the relevant difference in the distribution of the increments, which are no more Gaussian. 
