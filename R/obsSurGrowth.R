#' Estimate growth of individuals which survived between two censuses
#'
#' @param x dataframe of stem measurements
#' @param w column name of value from which to calculate growth in \code{x}
#' @param ind_id column name of individual IDs in code{x}
#' @param census_date column name of census dates in code{x}
#' @param census_date_1 column initial census of interval
#' @param census_date_2 column final census of interval
#'
#' @return total productivity from observed growth of survivors
#' 
#' @details These are stems which were observed as alive in both censuses
#' 
#' @export
#' 
obsSurGrowth <- function(x, w = "diam", ind_id = "stem_id", 
  census_date = "census_date", census_date_1, census_date_2) {

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

  # Find survivors
  si <- obsSur(x_fil, ind_id = ind_id, census_date = census_date, 
        census_date_1 = census_date_1, census_date_2 = census_date_2)

  # Filter to survivors
  x_si <- x_fil[x_fil[[ind_id]] %in% si, 
    c(ind_id, w, census_date)]

  # All stems have two measurements
  stopifnot(all(table(x_si[[ind_id]]) == 2))

  # Split into first and last census groups and order by stem ID
  x_cen1 <- x_si[x_si[[census_date]] == census_date_1,]
  x_cen2 <- x_si[x_si[[census_date]] == census_date_2,]
  x_cen1_w <- x_cen1[order(x_cen1[[ind_id]]), w]
  x_cen2_w <- x_cen2[order(x_cen2[[ind_id]]), w]

  # Calculate woody productivity, i.e. growth of recruit
  w_diff <- x_cen2_w - x_cen1_w

  # Sum of AGB growth per stem
  out <- sum(w_diff, na.rm = TRUE)

  # Return
  return(out)
}

