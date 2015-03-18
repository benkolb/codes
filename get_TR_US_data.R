
## Code by B. Kolb for Credit Spread News project, March 2015
# downloads US series and assembles them into .mat file

rm(list=ls())        # clear workspace

library(rdatamarket) # required for download from datamarket.com
library(xts)         # extended time series
library(pracma)      # for detrending

### Housekeeping
options(RCurlOptions = list(cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")))
# Set language to English (so "May2000" is recognised, instead of "Mai2000" etc.)
Sys.setlocale("LC_TIME", "English") 
# set working directory
if (Sys.info()[["nodename"]] == "SHAKAKHAN") setwd("C:/Users/Benedikt/Desktop/research/_projects/IIMPS_dsge/data")

### Download data and put into xts format
## Real GDP (Billions of Chained 2009 Dollars)
gdp_raw<- dmlist("https://datamarket.com/data/set/4hov/real-gross-domestic-product-per-capita#!ds=4hov")
gdp <- as.xts(gdp_raw$Value, as.yearqtr(gdp_raw$Quarter))
# consider growth rates
gdp <- setNames((gdp-lag(gdp))/lag(gdp)*400, "gdp")

## GDP Deflator (Gross Domestic Product - chain-type price index, % Change from Preceding Period)
def_raw <- dmlist("https://datamarket.com/data/set/4ssi/gross-domestic-product-chain-type-price-index#!ds=4ssi")
def <- as.xts(def_raw$Value, as.yearqtr(def_raw$Quarter))

## Federal funds rate (monthly)
ffr_raw <- dmlist("https://datamarket.com/data/set/1k1q/effective-federal-funds-rate-monthly-data#!ds=1k1q")
ffr <- as.xts(ffr_raw$Value, as.yearmon(ffr_raw$Month))
ffr <- apply.quarterly(ffr,'mean') # fill with previous value and average for quarters
index(ffr) <- as.yearqtr(index(ffr))


### Create and save dataset
dat <- na.omit(merge(gdp,def,ffr))

## Detrend GDP
# dat <- detrend(dat,'linear')

## Save as .csv
write.zoo(dat, "TR_US_data.csv", quote = FALSE, sep = ",")
