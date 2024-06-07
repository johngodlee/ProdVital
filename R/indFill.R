#' Create missing individual records
#'
#' @param x `r param_x()`
#' @param d vector of possible census dates for individuals in `x`
#' @param group `r param_group()`
#' @param census `r param_census()`
#'
#' @details 
#' `r details_group()`
#'
#' @examples
#' x <- data.frame(
#'   stem_id = "1",
#'   census_date = c(2000, 2002))
#' d = c(2000, 2001, 2002)
#' indFill(x, d, group = "stem_id", census = "census_date")
#'
#' @return dataframe of imputed individual records 
#' 
#' @export
#' 
indFill <- function(x, d, group, census) {

  # Convert potential tibble to dataframe
  x <- as.data.frame(x)

  # Generate unique individual IDs
  xint <- interaction(x[,group])

  # Split by group
  x_split <- split(x, xint, drop = TRUE)

  # For each individual:
  out <- do.call(rbind, lapply(x_split, function(y) {
    # Find first and last censuses for individual
    first_census <- min(y[[census]])
    last_census <- max(y[[census]])

    # Find missing censuses within first and last census
    missing_census <- d[
      !d %in% y[[census]] & 
        d < last_census & 
        d > first_census]

    # If any missing censuses:
    if (length(missing_census) > 0) {
      # Create empty rows and return
      new_meas <- unique(y[,group, drop = FALSE])
      new_meas <- new_meas[rep(1, length(missing_census)),, drop = FALSE]
      new_meas[[census]] <- missing_census
      rownames(new_meas) <- NULL
      new_meas
    }
  }))

  # Return
  return(out)
}


