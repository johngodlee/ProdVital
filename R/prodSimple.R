#' Calculate simple metrics of observed biomass change, across all
#'     census combinations for a single plot
#'
#' @param x dataframe of stem data from single plot
#' @param w column name of value from which to calculate growth in \code{x}
#' @param census_date column name of census dates in code{x}
#'
#' @return 
#' 
#' @examples
#' 
#' 
#' @export
#' 
prodSimple <- function(x, w = "diam", census_date = "census_date") {

  # Convert potential tibble to dataframe
  x <- as.data.frame(x)

  # Find all census dates
  census_date_all <- sort(unique(x[[census_date]]))

  # Stop if only one census
  if (length(census_date_all) == 1) { 
    stop("Plot has only one census")
  }

  # Create all pairwise combinations of censuses
  comb_list_pair <- combn(census_date_all, 2, simplify = FALSE)

  # For each pairwise combination of censuses
  out <- fastRbind(lapply(comb_list_pair, function(y) { 
    t0 <- y[1]
    tT <- y[2]
    int <- tT - t0
    x0 <- x[x[[census_date]] == t0,]
    xT <- x[x[[census_date]] == tT,]
    N0 <- nrow(x0)
    NT <- nrow(xT)
    dN <- NT - N0
    B0 <- sum(x0[[w]], na.rm = TRUE)
    BT <- sum(xT[[w]], na.rm = TRUE)
    dB <- BT - B0
    dB_rel <- dB / B0
    data.frame(
      t0, tT, int,
      N0, NT, dN, 
      B0, BT, dB, dB_rel)
  }))

  # Return
  return(out)
}







