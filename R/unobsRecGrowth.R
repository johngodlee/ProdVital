#' Estimate growth of unobserved recruits between two censuses
#'
#' @param x dataframe of stem measurements
#' @param w column name of value from which to calculate growth in \code{x}
#' @param ind_id column name of individual IDs in code{x}
#' @param diam column name of stem diameters in code{x}
#' @param census_date column name of census dates in code{x}
#' @param min_size_class vector of length two containing range of minimum 
#'     diameter size class used to estimate median growth rate
#' @param census_date_1 column initial census of interval
#' @param census_date_2 column final census of interval
#' @param w_min the minimum value of \code{w} expected to be encountered in 
#'     the plot given the minimum diameter threshold. Default assumes minimum 
#'     diameter threshold
#'
#' @return total productivity from unobserved growth of recruits which died
#' 
#' @details These are stems which recruit in sometime after the previous 
#'     census, but die before the next census.
#' 
#' @export
#' 
unobsRecGrowth <- function(x, w = "diam", ind_id = "stem_id", diam = "diam",
  census_date = "census_date", min_size_class = c(5, 10), 
  census_date_1, census_date_2, w_min = min_size_class[1]) {

  # Convert potential tibble to dataframe
  x <- as.data.frame(x)

  # Find census interval
  int <- as.numeric(census_date_2) - as.numeric(census_date_1)

  # Stop is interval is negative
  if (int < 0) { 
    stop("census dates must be in chronological order")
  }

  # Subset to censuses of interest
  x_fil <- x[x[[census_date]] %in% c(census_date_1, census_date_2),]

  # Estimate number of unobserved recruits in census interval
  unobs_rec <- unobsRec(x_fil, census_date = census_date, 
    census_date_1 = census_date_1, census_date_2 = census_date_2,
    return_interm = FALSE)

  # Find observed survivors in final census
  si <- obsSur(x_fil, ind_id = ind_id, census_date = census_date, 
    census_date_1 = census_date_1, census_date_2 = census_date_2)
  x_si <- x_fil[x_fil[[ind_id]] %in% si, 
    c(ind_id, w, census_date, diam)]

  # Subset individuals to smallest size class
  x_si_min <- x_si[x_si[[diam]] >= min_size_class[1] & 
    x_si[[diam]] < min_size_class[2] & 
    !is.na(x_si[[diam]]),]

  # Find median growth rate of survivors in smallest diam size class 
  # Between two censuses
  si_growth_median <- median(
    indGrowth(x_si_min, w = w, ind_id = ind_id, 
      census_date_1 = census_date_1, census_date_2 = census_date_2))
  
  # Calculate woody productivity, i.e. change of recruit
  w_diff <- (w_min + (w_min * si_growth_median)) - w_min

  # Multiply by 1/3 (recruit at 1/3 interval, die at 2/3 interval)
  # Multiply by number of unobserved recruits
  out <- w_diff * 1/3 * unobs_rec

  # Return
  return(out)
}

