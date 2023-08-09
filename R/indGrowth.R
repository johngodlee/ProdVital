#' Calculate growth increments 
#' 
#' `r descrip_table()` to calculate the growth increment of each individual
#' given a measure of individual size (`w`).
#'
#' @param x `r param_x()`
#' @param t0 `r param_t0()`
#' @param tT `r param_tT()`
#' @param w `r param_w()`
#' @param group `r param_group()`
#' @param census `r param_census()`
#' 
#' @return 
#' Named vector of individual growth increment. Names are the concatenation 
#' (sep = `.`) of all grouping columns provided in `group`.
#' 
#' @details 
#' `r details_group()`
#' 
#' `r details_t0_tT()`
#'
#' @examples
#' data(bicuar)
#' 
#' indGrowth(bicuar, t0 = "2019", tT = "2021", w = "agb", 
#'   group = "stem_id", census = "census_date")
#' 
#' @export
#' 
indGrowth <- function(x, t0, tT, w, group, census) {

  # Assert as dataframe
  x <- as.data.frame(x)

  # Calculate census interval
  int <- as.numeric(tT) - as.numeric(t0)

  # Stop is interval is negative
  if (int < 0) { 
    stop("Census dates must be in chronological order")
  }

  # Stop if any columns not recognised
  if (any(!c(w, group, census) %in% names(x))) {
    stop("Some columns not present in x")
  }

  # Subset to censuses of interest
  x_fil <- x[x[[census]] %in% c(t0, tT),]

  # Stop if a census is missing
 if (!all(c(t0, tT) %in% x[[census]])) {
    stop("One or more censuses is missing from data")
  }

  # Find individuals present in both censuses
  si <- obsID(x_fil, type = "sur", group = group, census = census,
    t0 = t0, tT = tT)

  # Filter to survivors
  x_si <- merge(x_fil, si) 

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

