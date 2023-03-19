#' Calculate the observed loss to mortality between two censuses 
#'     for a single plot
#'
#' @param x dataframe of SEOSAW stem data from single plot
#' @param ind_id column name of ind IDs 
#' @param census_date column name of census dates 
#' @param census_date_1 column initial census of interval
#' @param census_date_2 column final census of interval
#' 
#' @return observed loss to mortality
#'  
#' @export
#' 
obsMorLoss <- function(x, w = "diam", ind_id = "stem_id", 
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

  # Find individuals which died between two censuses
  di <- obsMor(x_fil, ind_id = ind_id, census_date = census_date,
    census_date_1 = census_date_1, census_date_2 = census_date_2)
  x_di <- x_fil[x_fil[[ind_id]] %in% di & 
      x_fil[[census_date]] == census_date_1, 
    c(ind_id, census_date, w)]

  # Calculate total loss
  out <- sum(x_di[[w]], na.rm = TRUE)

  # Should never be negative
  if (out < 0) { stop("No observed mortality should be <0") }

  # Return
  return(out)
}

