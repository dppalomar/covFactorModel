#' @title Sector Information of Stocks
#' 
#' @description Obtains the sector information of a set of stocks.
#'  
#' @param stocks vector containing the stock names
#' @param stock_sector_database matrix containing the pairing of stocks and sectors (first column contains the stock names and second column sector names)
#' @return A list with the following components:
#'         \item{\code{stock_sector_info}}{integer vector with each number being an index representing one specific sector}
#'         \item{\code{sectors}}{string vector containing sector names sorted corresponding to their sector index}
#' @author ZHOU Rui and Daniel P. Palomar
#' @examples
#' mystocks <- c("AAPL",  "ABBV", "AET", "AMD", "APD", "AA","AVY", "A", "ADI", "IBM")
#' getSectorInfo(mystocks)
#' @export getSectorInfo
getSectorInfo <- function(stocks, stock_sector_database) {
  # load default database if needed
  if (missing(stock_sector_database))
    data(stock_sector_database, envir = environment())
  
  # match sectors
  mask <- match(stocks, stock_sector_database[, 1])
  if (anyNA(mask)) {
    stocks_no_sector <- stocks[is.na(mask)]
    stop("Cannot find sector info to these stock symbols: ", paste(stocks_no_sector, collapse = ", "))
  }
  
  # extract sector info
  sectors_full <- stock_sector_database[mask, 2]
  sectors <- unique(sectors_full)
  stock_sector_info <- match(sectors_full, sectors)
  
  # assign name
  names(stock_sector_info) <- stocks
  names(sectors) <- c(1:length(sectors))
  
  return(list(
    "stock_sector_info" = stock_sector_info,
    "sectors" = sectors
  ))
}  