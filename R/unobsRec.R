#' Estimate unobserved number of recruits according to Talbot et al. (2014)
#' 
#' `r descrip_table()` to estimate the number of recruits which recruited and
#' died within a single census interval.
#'
#' @param x `r param_x()`
#' @param t0 `r param_t0()`
#' @param tT `r param_tT()`
#' @param group `r param_group()`
#' @param census `r param_census()`
#'
#' @return 
#' `details_obs_sum(un = TRUE)` unobserved recruits.
#' 
#' @details
#' `r details_group()`
#' 
#' @examples
#' data(bicuar)
#' 
#' unobsRec(bicuar, "2019", "2021", group = "stem_id", census = "census_date")
#' 
#' @export
#' 
unobsRec <- function(x, t0, tT, group, census) { 

  # Convert potential tibble to dataframe
  x <- as.data.frame(x)

  # Number of living stems in plot at t2
  N <- sum(x[[census]] == tT)

  # Number of living stems in plot at t1
  n0 <- sum(x[[census]] == t0)

  # Number of deaths between censuses
  nt_d <- nrow(obsID(x, t0 = t0, tT = tT, type = "mor", 
      group = group, census = census))

  # Number of recruits between censuses
  nt_r <- nrow(obsID(x, t0 = t0, tT = tT, type = "rec", 
      group = group, census = census))

  # Census interval
  int <- as.numeric(tT) - as.numeric(t0)

  # Mean annual mortality rate in plot
  M <- nt_d / n0 / int

  # Mean annual recruitment rate in plot
  R <- nt_r / n0 / int

  # Estimate unobserved recruitss
  Ur <- N*M*R*int

  # Return
  return(Ur)
}

