#' Calculate instantaneous turnover rates - Kohyama et al. (2019)
#'
#' Used in \code{kohyamaProd()}.
#' 
#' @param y initial value
#' @param z final value
#' @param t census interval
#'
#' @keywords internal
#' @noRd
#' 
turnoverEst <- function(y, z, t) {
  f <- function(rho) {
    sum(y * exp(-rho * t) - z)
  }
  df <- function(rho) {
    sum(-t * y * exp(-rho * t))
  }
  # Newton-Rapton iteration
  rho <- 0.02
  precision <- 1.0e-12 # to stop iteration
  change <- precision + 1.0

  while (change > precision) {
    rho2 <- rho - f(rho) / df(rho)
    change <- abs(rho2 - rho)
    rho <- rho2
  }
  return(rho)
}

