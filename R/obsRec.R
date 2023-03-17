#' Find IDs of recruits between two censuses for a single plot
#'
#' @param x dataframe of SEOSAW stem data from single plot
#' @param ind_id column name of ind IDs 
#' @param census_date column name of census dates 
#' @param census_date_1 column initial census of interval
#' @param census_date_2 column final census of interval
#'
#' @return vector of individual ID values for stems which 
#'     recruited between censuses 
#' 
#' @export
#' 
obsRec <- function(x, ind_id = "stem_id", census_date = "census_date", 
  census_date_1, census_date_2) {

  # Convert potential tibble to dataframe
  x <- as.data.frame(x)

  # Subset to initial census
  x0 <- x[x[[census_date]] == census_date_1, ind_id]

  # Subset to final census
  xt <- x[x[[census_date]] == census_date_2, ind_id]

  # Find stem IDs in xt but not in x0
  out <- xt[!xt %in% x0]

  # Return
  return(out)
}

