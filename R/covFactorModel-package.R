#'covFactorModel: Covariance matrix estimation via factor models
#'
#'Estimation of covariance matrix via factor models with application
#'to financial data. Factor models decompose the asset returns into an 
#'exposure term to some factors and a residual idiosynchratic component. The 
#'resulting covariance matrix contains a low-rank term corresponding to the 
#'factors and another full-rank term corresponding to the residual component. 
#'This package provides a function to separate the data into the factor 
#'component and residual component, as well as to estimate the corresponding 
#'covariance matrix. Different kind of factor models are considered, namely, 
#'macroeconomic factor models and statistical factor models. The estimation
#'of the covariance matrix accepts different kinds of structure on the 
#'residual term: diagonal structure (implying that residual component is 
#'uncorrelated) and block diagonal structure (allowing correlation within 
#'sectors). The package includes a built-in database containing stock symbol
#'and their sectors.
#'
#' @section Functions:
#' \code{\link{factorModel}}, \code{\link{covFactorModel}}, \code{\link{getSectorInfo}}
#'
#' @section Help:
#' For a quick help see the README:
#' \href{https://rawgit.com/dppalomar/covFactorModel/master/README.html}{GitHub-README}.
#'  
#' For more details see the vignette:
#' \href{https://rawgit.com/dppalomar/covFactorModel/master/vignettes/covFactorModel-vignette.html}{GitHub-html-vignette}, and
#' \href{https://rawgit.com/dppalomar/covFactorModel/master/vignettes/covFactorModel-vignette.pdf}{GitHub-pdf-vignette}.
#' 
#' @author Rui ZHOU and Daniel P. Palomar
#'
#' @references
#' R. S. Tsay, \emph{Analysis of Financial Time Series.} John Wiley & Sons, 2005.
#' 
#' @docType package
#' @name covFactorModel-package
NULL
