#' Calculate loss from individuals which died 
#' 
#' `r descrip_table()` `r descrip_gro("calculate", "loss", "died")
#'
#' @param x `r param_x()`
#' @param t0 `r param_t0()`
#' @param tT `r param_tT()`
#' @param w `r param_w()`
#' @param group `r param_group()`
#' @param census `r param_census()`
#' 
#' @return 
#' `r details_obs_sum()` loss from individuals which died.
#' 
#' @details
#' `r details_group()`
#'
#' @examples
#' data(bicuar)
#'
#' obsMorLoss(bicuar, "2019", "2021", w = "agb", 
#'   group = "stem_id", census = "census_date")
#'  
#' @export
#' 
obsMorLoss <- function(x, t0, tT, w, group, census) {

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

  # Find individuals which died between two censuses
  di <- obsID(x_fil, type = "mor", group = group, census = census,
    t0 = t0, tT = tT)

  # Filter to deaths
  x_di <- merge(x_fil, di)
  x_loss <- x_di[x_di[[census]] == t0, w]

  # Calculate total loss
  out <- sum(x_loss, na.rm = TRUE)

  # Should never be negative
  if (out < 0) { stop("No observed mortality should be <0") }

  # Return
  return(out)
}

