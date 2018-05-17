#### This file is for examples

### Example 1 -- build factor model with Index
### choose 10 stocks from S&P500, and use S&P500 Index as explicit factor 
library(xts)
library(quantmod)

## set begin-end date and stock namelist
begin_date <- "2017-01-01"
end_date <- "2017-06-30"
stock_namelist <- c("A", "APD", "AAL", "AAP", "AAPL", "ADM", "AEE", "AFL", "APA", "AVB")

## download data from internet
dataSet <- xts()
for (stock_index in 1:length(stock_namelist)) {
  dataSet <- merge(dataSet, getSymbols(stock_namelist[stock_index], from = begin_date, to = end_date, auto.assign = FALSE)[, 4])
}
SP500_index <- getSymbols("^GSPC", from = begin_date, to = end_date, auto.assign = FALSE)[, 4]
indexClass(dataSet) <- "Date"
indexClass(SP500_index) <- "Date"

## rename data columns
colnames(dataSet) <- stock_namelist
colnames(SP500_index) <- "INDEX"

## use SP500 index as explicit factor
X <- diff(log(dataSet), lag = 1, na.pad = FALSE)
factor_index <- diff(log(SP500_index), lag = 1, na.pad = FALSE)
factor_model <- removeFactors(X = X, type = "E", explicit_factor = factor_index)



### Example 2 -- use S&P500 index to analyze fund performance (assuming rf = 0)
### choose 5 ETF which are:
### SPY - SPDR S&P 500 ETF (index tracking)
### XIVH - Velocity VIX Short Vol Hedged (high alpha)
### SPHB - PowerShares S&P 500 High Beta Portfolio (high Beta)
### SPLV - PowerShares S&P 500 Low Volatility Portfolio (low Beta)
### USMV - iShares Edge MSCI Min Vol USA ETF 
### JKD - iShares Morningstar Large-Cap ETF 
library(xts)
library(quantmod)

## set begin-end date and stock namelist
begin_date <- "2016-10-01"
end_date <- "2017-06-30"
stock_namelist <- c("SPY","XIVH", "SPHB", "SPLV", "USMV", "JKD")

## download data from internet
dataSet <- xts()
for (stock_index in 1:length(stock_namelist)) {
  dataSet <- merge(dataSet, getSymbols(stock_namelist[stock_index], from = begin_date, to = end_date, auto.assign = FALSE)[, 4])
}
SP500_index <- getSymbols("^GSPC", from = begin_date, to = end_date, auto.assign = FALSE)[, 4]
indexClass(dataSet) <- "Date"
indexClass(SP500_index) <- "Date"

## rename data columns
colnames(dataSet) <- stock_namelist
colnames(SP500_index) <- "INDEX"

## use SP500 index as explicit factor
X <- diff(log(dataSet), lag = 1, na.pad = FALSE)
factor_index <- diff(log(SP500_index), lag = 1, na.pad = FALSE)
factor_model <- removeFactors(X = X, type = "E", explicit_factor = factor_index)

## show Alpha and Beta result in plot
result <- rbind(t(factor_model$const)*1e3, t(factor_model$loading))
rownames(result) <- c("Alpha*1000", "Beta")
barplot(height = result,
        col = c("red", "blue"),
        beside = TRUE,
        legend.text = TRUE)


### Example 3 -- Fama and French Model
### choose 10 stocks from S&P500, and use Fama-French Factors 
library(xts)
library(quantmod)

## set begin-end date and stock namelist
begin_date <- "2017-01-01"
end_date <- "2017-06-30"
stock_namelist <- c("A", "APD", "AAL", "AAP", "AAPL", "ADM", "AEE", "AFL", "APA", "AVB")

## download data from internet
# download stocks price
dataSet <- xts()
for (stock_index in 1:length(stock_namelist)) {
  dataSet <- merge(dataSet, getSymbols(stock_namelist[stock_index], from = begin_date, to = end_date, auto.assign = FALSE)[, 4])
}
indexClass(dataSet) <- "Date"
# download Fama- French factors
url <- "http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Research_Data_Factors_daily_CSV.zip"
temp <- tempfile()
download.file(url, temp, method = "libcurl", mode = "wb")
unzip(temp, "F-F_Research_Data_Factors_daily.CSV")
mydata <- read.csv("F-F_Research_Data_Factors_daily.CSV", skip = 4)
unlink(temp)
mydata <- mydata[-dim(mydata)[1], ]
fama_lib <- xts(x = mydata[, -1], order.by = as.Date(paste(mydata[, 1]), "%Y%m%d"))

## rename data columns
colnames(dataSet) <- stock_namelist

## Fama-French Model
X <- diff(log(dataSet), lag = 1, na.pad = FALSE)
X <- X - as.matrix(fama_lib[index(X), 4]/100) %*% matrix(data = 1, nrow = 1, ncol = length(stock_namelist))
factors_fama <- fama_lib[index(X), -4]/100
factor_model <- removeFactors(X = X, type = "E", explicit_factor = factors_fama)


### Example 4 -- clustering VS. prior sector info
### cluster stocks with 1-factor model residual correlation matrix and compare it with prior sector information
library(xts)
library(quantmod)
library(pheatmap)

## set begin-end date and stock namelist
begin_date <- "2016-01-01"
end_date <- "2017-06-30"
stock_namelist <- c("AAPL",  "ABBV", "AET", "AMD", "APD", "AA","CF", "A", "ADI", "IBM")

## load and achieve prior sector info
data(sectors_info)
prior_sector_info <- getSectorInfo(stock_list = stock_namelist, sector_info = sectors_info)

## download data from internet
dataSet <- xts()
for (stock_index in 1:length(stock_namelist)) {
  dataSet <- merge(dataSet, getSymbols(stock_namelist[stock_index], from = begin_date, to = end_date, auto.assign = FALSE)[, 4])
}
SP500_index <- getSymbols("^GSPC", from = begin_date, to = end_date, auto.assign = FALSE)[, 4]
indexClass(dataSet) <- "Date"
indexClass(SP500_index) <- "Date"

## rename data columns and add sector info
colnames(dataSet) <- paste(stock_namelist, '_', paste(prior_sector_info$stock_sector_info), sep = '')
colnames(SP500_index) <- "INDEX"

## use SP500 index as explicit factor
X <- diff(log(dataSet), lag = 1, na.pad = FALSE)
factor_index <- diff(log(SP500_index), lag = 1, na.pad = FALSE)
factor_model <- removeFactors(X = X, type = "E", explicit_factor = factor_index)

## use 1-factor model residual correlation matrix to cluster stocks
cor_residual <- cor(factor_model$perturbations)
pheatmap(mat = cor_residual)
