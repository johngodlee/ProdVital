#' Estimate unobserved number of recruits according to Talbot et al. (2014)
#'
#' @param x dataframe of SEOSAW stem data from single plot
#' @param ind_id column name of individual IDs 
#' @param census_date column name of census dates 
#' @param census_date_1 column initial census of interval
#' @param census_date_2 column final census of interval
#' @param return_interm logical, return intermediate parameters as well?
#'
#' @return estimate of unobserved recruitment, 
#'     optionally as a dataframe with intermediate parameters.
#' 
#' @export
#' 
unobsRec <- function(x, ind_id = "stem_id", census_date = "census_date", 
  census_date_1, census_date_2, return_interm = FALSE) { 

  # Convert potential tibble to dataframe
  x <- as.data.frame(x)

  # Number of living stems in plot at t2
  N <- sum(x[[census_date]] == census_date_2)

  # Number of living stems in plot at t1
  n0 <- sum(x[[census_date]] == census_date_1)

  # Number of deaths between censuses
  nt_d <- length(obsMor(x, ind_id = ind_id, census_date = census_date, 
    census_date_1 = census_date_1, census_date_2 = census_date_2))

  # Number of recruits between censuses
  nt_r <- length(obsRec(x, ind_id = ind_id, census_date = census_date, 
    census_date_1 = census_date_1, census_date_2 = census_date_2))

  # Census interval
  int <- as.numeric(census_date_2) - as.numeric(census_date_1)

  # Mean annual mortality rate in plot
  M <- nt_d / n0 / int

  # Mean annual recruitment rate in plot
  R <- nt_r / n0 / int

  # Estimate unobserved recruitss
  Ur <- N*M*R*int

  # Optionally return dataframe with intermediate parameters
  if (return_interm) { 
    data.frame(N, n0, nt_d, nt_r, int, M, R, Ur)
  } else { 
    return(Ur)
  }
}

