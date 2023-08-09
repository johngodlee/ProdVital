#' Estimate growth of unobserved recruits between two censuses
#'
#' `r descrip_table()` `r descrip_gro("estimate", "growth", "recruited and died")`
#'
#' @param x `r param_x()`
#' @param t0 `r param_t0()`
#' @param tT `r param_tT()`
#' @param w `r param_w()`
#' @param group `r param_group()`
#' @param census `r param_census()`
#' @param diam `r param_diam()`
#' @param min_size_class `r param_min_size_class()`
#' @param w_min_diam `r param_w_min_diam()`
#'
#' @return 
#' `r details_obs_sum(un = TRUE)` growth from recruits which died.
#' 
#' @details 
#' `r details_group()`
#' 
#' @examples
#' data(bicuar)
#' 
#' unobsRecGrowth(bicuar, "2019", "2021", w = "diam", group = "stem_id", 
#'   census = "census_date", diam = "diam", min_size_class = c(5, 10), 
#'   w_min_diam = 5)
#' 
#' @importFrom stats median
#'
#' @export
#' 
unobsRecGrowth <- function(x, t0, tT, w, group, census, diam, 
  min_size_class, w_min_diam) {

  # Stop if minimum size class is malformed
  if (length(min_size_class) != 2 | !is.numeric(min_size_class)) {
    stop("min_size_class must be a numeric vector of length two")
  }

  # Convert potential tibble to dataframe
  x <- as.data.frame(x)

  # Find census interval
  int <- as.numeric(tT) - as.numeric(t0)

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
  if (any(!c(w_min_diam, w, group, diam, census) %in% names(x))) {
    stop("Some columns not present in x")
  }

  # Subset to censuses of interest
  x_fil <- x[x[[census]] %in% c(t0, tT),]

  # Estimate number of unobserved recruits in census interval
  unobs_rec <- unobsRec(x_fil, t0 = t0, tT = tT, 
    group = group, census = census)

  # Find observed survivors in final census
  si <- obsID(x_fil, type = "sur", group = group, census = census, 
    t0 = t0, tT = tT)
  x_si <- merge(x_fil, si)

  # Subset individuals to smallest size class
  x_si_min <- x_si[x_si[[diam]] >= min_size_class[1] & 
    x_si[[diam]] < min_size_class[2] & 
    !is.na(x_si[[diam]]),]

  # Find median growth rate of survivors in smallest diam size class 
  # Between two censuses
  si_growth_median <- stats::median(
    indGrowth(x_si_min, w = w, group = group, census = census,
      t0 = t0, tT = tT), 
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


