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
#' data(bicuar_clean)
#' 
#' obsSurGrowth(bicuar_clean, "2019", "2021", w = "diam", 
#'   group = c("plot_id", "stem_id"), census = "census_date")
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
  stopifnot(all(table(as.character(interaction(x_si[,group, drop = TRUE]))) == 2))

  # Filter to first and last census
  x_si0 <- x_si[x_si[[census]] == t0,]
  x_sit <- x_si[x_si[[census]] == tT,]

  # order by group
  x_si0_ord <- x_si0[order(interaction(x_si0[,group])),]
  x_sit_ord <- x_sit[order(interaction(x_sit[,group])),]

  # diff
  out <- x_sit_ord[[w]] - x_si0_ord[[w]]

  # Add names
  names(out) <- interaction(x_si0_ord[,group])

  # Return
  return(out)
}

