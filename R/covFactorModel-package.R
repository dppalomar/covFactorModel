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
#' \href{https://rawgit.com/dppalomar/covFactorModel/master/README.html}{GitHub-README} and
#' \href{https://cran.r-project.org/web/packages/covFactorModel/README.html}{CRAN-README}.
#'  
#' For more details see the vignette:
#' \href{https://rawgit.com/dppalomar/covFactorModel/master/vignettes/covFactorModel-vignette.html}{GitHub-html-vignette},
#' \href{https://rawgit.com/dppalomar/covFactorModel/master/vignettes/covFactorModel-vignette.pdf}{GitHub-pdf-vignette}, and
#' \href{https://cran.r-project.org/web/packages/covFactorModel/vignettes/covFactorModel-vignette.pdf}{CRAN-pdf-vignette}.
#' 
#' @author Rui ZHOU and Daniel P. Palomar
#'
#' @references
#' R. S. Tsay, \emph{Analysis of Financial Time Series.} John Wiley & Sons, 2005.
#' 
#' @docType package
#' @name covFactorModel-package
NULL
