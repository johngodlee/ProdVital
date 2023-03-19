#' Estimate growth of unob
#' Estimate growth of unobserved recruits between two censuses
#'
#' @param x dataframe of stem measurements
#' @param w column name of value from which to calculate growth in \code{x}
#' @param ind_id column name of individual IDs in code{x}
#' @param diam column name of stem diameters in code{x}
#' @param census_date column name of census dates in code{x}
#' @param min_size_class vector of length two containing range of minimum 
#'     diameter size class used to estimate median growth rate
#' @param w_min_diam vector of length 1, column name of estimated value of 
#'     \code{w} in \code{x} at the minimum diameter threshold of the plot. 
#'     Defaults to lower value of minimum diameter size class
#' @param census_date_1 column initial census of interval
#' @param census_date_2 column final census of interval
#'
#' @return total productivity from unobserved growth of recruits which died
#' 
#' @details These are stems which recruit in sometime after the previous 
#'     census, but die before the next census.
#' 
#' @export
#' 
unobsRecGrowth <- function(x, w = "diam", ind_id = "stem_id", 
  diam = "diam", census_date = "census_date", 
  min_size_class = c(5, 10), w_min_diam = min_size_class[1], 
  census_date_1, census_date_2) {

  # Stop if minimum size class is malformed
  if (length(min_size_class) != 2 | !is.numeric(min_size_class)) {
    stop("min_size_class must be a numeric vector of length two")
  }

  # Convert potential tibble to dataframe
  x <- as.data.frame(x)

  # Find census interval
  int <- as.numeric(census_date_2) - as.numeric(census_date_1)

  # Stop is interval is negative
  if (int < 0) { 
    stop("census dates must be in chronological order")
  }

  # If w_min_diam is numeric, add as column
  if (is.numeric(w_min_diam)) {
    if (length(w_min_diam) > 1) {
      stop("Numeric w_min_diam must be of length one")
    }
    x$w_min_diam <- w_min_diam
    w_min_diam <- "w_min_diam"
  } 

  # Stop if any columns not recognised
  if (any(!c(w_min_diam, w, ind_id, diam, census_date) %in% names(x))) {
    stop("Some columns not present in x")
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
    c(ind_id, w, census_date, diam, w_min_diam)]

  # Subset individuals to smallest size class
  x_si_min <- x_si[x_si[[diam]] >= min_size_class[1] & 
    x_si[[diam]] < min_size_class[2] & 
    !is.na(x_si[[diam]]),]

  # Find median growth rate of survivors in smallest diam size class 
  # Between two censuses
  si_growth_median <- median(
    indGrowth(x_si_min, w = w, ind_id = ind_id, 
      census_date_1 = census_date_1, census_date_2 = census_date_2), 
    na.rm = TRUE)

  # Find median w of those stems at minimum diameter threshold
  si_w_min_diam_median <- median(x_si_min[[w_min_diam]])
  
  # Calculate median woody productivity, i.e. change of recruit
  w_diff <- si_w_min_diam_median * si_growth_median

  # Multiply by 1/3 (recruit at 1/3 interval, die at 2/3 interval)
  # Multiply by number of unobserved recruits
  out <- w_diff * 1/3 * unobs_rec

  # If missing, assume zero
  if (length(out) == 0) {
    out <- 0
  }

  # Return
  return(out)
}


