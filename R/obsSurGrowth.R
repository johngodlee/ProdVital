#' Calculate growth from individuals which survived between two censuses
#'
#' `r descrip_table()` `r descrip_gro("calculate", "growth", "survived")`
#'
#' @param x `r param_x()`
#' @param t0 `r param_t0()`
#' @param tT `r param_tT()`
#' @param w `r param_w()`
#' @param group `r param_group()`
#' @param census `r param_census()`
#'
#' @return 
#' `r details_obs_sum()` growth from individuals which survived.
#' 
#' @details
#' `r details_group()`
#'
#' @examples
#' data(bicuar)
#' 
#' obsSurGrowth(bicuar, "2019", "2021", w = "diam", 
#'   group = "stem_id", census = "census_date")
#' 
#' @export
#' 
obsSurGrowth <- function(x, t0, tT, w, group, census) {

  # Convert potential tibble to dataframe
  x <- as.data.frame(x)

  # Find census interval
  int <- as.numeric(tT) - as.numeric(t0)

  # Stop is interval is negative
  if (int < 0) { 
    stop("census dates must be in chronological order")
  }

  # Subset to censuses of interest
  x_fil <- x[x[[census]] %in% c(t0, tT),]

  # Find survivors
  si <- obsID(x_fil, type = "sur", group = group, census = census, 
        t0 = t0, tT = tT)

  # Filter to survivors
  x_si <- merge(x_fil, si)

  # All stems have two measurements
  stopifnot(all(table(interaction(x_si[,group, drop = FALSE])) == 2))

  # Split into first and last census groups and order by stem ID
  x_cen1 <- x_si[x_si[[census]] == t0,]
  x_cen2 <- x_si[x_si[[census]] == tT,]
  x_cen1_w <- x_cen1[order(interaction(x_cen1[, group, drop = FALSE])), w]
  x_cen2_w <- x_cen2[order(interaction(x_cen2[, group, drop = FALSE])), w]

  # Calculate woody productivity, i.e. growth of stems
  w_diff <- x_cen2_w - x_cen1_w

  # Sum of AGB growth per stem
  out <- sum(w_diff, na.rm = TRUE)

  # Return
  return(out)
}

