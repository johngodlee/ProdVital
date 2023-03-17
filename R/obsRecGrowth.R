#' Estimate growth of individuals which recruited between two censues
#'
#' @param x dataframe of stem measurements
#' @param w column name of value from which to calculate growth in \code{x}
#' @param ind_id column name of individual IDs in code{x}
#' @param diam column name of stem diameters in code{x}
#' @param census_date column name of census dates in code{x}
#' @param growth_percentile percentile of growth rate used to estimate growth 
#'     rate of recruits for backward extrapolation
#' @param min_size_class vector of length two containing range of minimum 
#'     diameter size class used to estimate median growth rate of recruits
#' @param rec_method character string describing method used to estimate 
#'     growth of recruits
#' @param min_diam_thresh minimum diameter threshold of plot
#' @param census_date_1 column initial census of interval
#' @param census_date_2 column final census of interval
#' @param w_min the minimum value of \code{w} expected to be encountered in 
#'     the plot given the minimum diameter threshold. Default assumes minimum 
#'     diameter threshold
#'
#' @return total productivity from observed growth of recruits
#' 
#' @details These are stems which were observed as a recruit in the next 
#'     census.
#' @export
#' 
obsRecGrowth <- function(x, w = "diam", ind_id = "stem_id", diam = "diam",
  census_date = "census_date",  growth_percentile = 0.86, 
  min_size_class = c(5, 10),
  rec_method = c("zero", "min_diam_thresh", "extrap"), 
  min_diam_thresh, census_date_1, census_date_2, w_min = min_size_class[1]) { 

  # Stop if methods not recognised
  if (any(!rec_method %in% c("zero", "min_diam_thresh", "extrap"))) {
    stop("Invalid recruit growth estimation method. Methods must be either 'zero', 'min_diam_thresh', 'extrap'")
  }

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

  # Find recruits
  ri <- obsRec(x_fil, ind_id = ind_id, census_date = census_date,
      census_date_1 = census_date_1, census_date_2 = census_date_2)

  # Create list to hold output
  method_list <- vector("list", length(rec_method))
  names(method_list) <- rec_method

  # If there are recruits:
  if (length(ri) > 0) {
    # Filter to recruits
    x_rec <- x_fil[x_fil[[ind_id]] %in% ri,]

    # Method: Assume grew from 0 since previous census (zero)
    if ("zero" %in% rec_method) {
      # Filter recruits to final census
      x_rec_fil <- x_rec[x_rec[[census_date]] == census_date_2,]

      # Calculate sum of growth
      BD <- sum(x_rec_fil[[w]], na.rm = TRUE)

      # Create list
      method_list[["zero"]] <- BD
    }

    # Method: Assume grew from min. diam. thresh. (min_diam_thresh)
    if ("min_diam_thresh" %in% rec_method) {

      # Filter recruits to final census
      x_rec_fil <- x_rec[x_rec[[census_date]] == census_date_2,]

      # For each recruit make w the min. diam. thresh equivalent
      x_rec_fil$w_old <- w_min 

      # Find change
      BD <- sum(x_rec_fil[[w]] - x_rec_fil$w_old, na.rm = TRUE)

      # Create list
      method_list[["min_diam_thresh"]] <- BD
    }

    # Extrap. diam. to first census using percentile of small stem growth rate (extrap)
    if ("extrap" %in% rec_method) {

      # Filter recruits to final census
      x_rec_fil <- x_rec[x_rec[[census_date]] == census_date_2,]

      # Find survivors
      si <- obsSur(x_fil, ind_id = ind_id, census_date = census_date, 
        census_date_1 = census_date_1, census_date_2 = census_date_2)

      # Filter to survivors
      x_si <- x_fil[x_fil[[ind_id]] %in% si, c(ind_id, w, diam, census_date)]

      # Find stems in smallest diameter size class
      x_si_min <- x_si[x_si[[diam]] >= min_size_class[1] & 
        x_si[[diam]] < min_size_class[2],]

      # Find percentile growth rate of survivors in smallest diam size class 
      si_diam_growth_median <- quantile(
        indGrowth(x_si_min, w = diam, ind_id = ind_id, 
          census_date_1 = census_date_1, census_date_2 = census_date_2), 
        growth_percentile)

      # For each stem, extrapolate back the stem diameter 
      x_rec_fil$diam_cen1 <- x_rec_fil$diam - si_diam_growth_median

      # Filter to extrapolated diameters > min. diam. thresh
      x_rec_gemin <- x_rec_fil[x_rec_fil$diam_cen1 > min_diam_thresh,]

      # Find change
      BD <- sum(x_rec_gemin[[w]], na.rm = TRUE)

      # Create list
      method_list[["extrap"]] <- BD
    }
  }

  # If only one method, don't return nested list
  if (length(rec_method) == 1) {
    method_list <- method_list[[1]]
  }

  # Return
  return(method_list)
}

