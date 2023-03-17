#' Calculate stem growth increment for a single group over a single census interval
#'
#' @param x dataframe of SEOSAW stem data 
#' @param w column name of value from which to calculate growth in \code{x}
#' @param ind_id column name of individual IDs 
#' @param census_date column name of census dates 
#' @param census_date_1 column initial census of interval
#' @param census_date_2 column final census of interval
#'
#' @return vector of individual growth increment
#' 
#' @export
#' 
indGrowth <- function(x, w = "diam", ind_id = "stem_id",
  census_date = "census_date", census_date_1, census_date_2) {

  # Convert to dataframe
  x <- as.data.frame(x)

  # Calculate census interval
  int <- as.numeric(census_date_2) - as.numeric(census_date_1)

  # Stop is interval is negative
  if (int < 0) { 
    stop("census dates must be in chronological order")
  }

  # Subset to censuses of interest
  x_fil <- x[x[[census_date]] %in% c(census_date_1, census_date_2),]

  # Find individuals present in both censuses
  si <- obsSur(x_fil, ind_id = ind_id, census_date = census_date,
    census_date_1 = census_date_1, census_date_2 = census_date_2)

  # Filter to survivors
  x_si <- x_fil[x_fil[[ind_id]] %in% si,]

  # Split by stem ID
  x_split <- split(x_si, x_si[[ind_id]])

  # All should have two rows
  stopifnot(all(unlist(lapply(x_split, nrow)) == 2))

  # Calculate growth rate of each stem
  out <- unlist(lapply(x_split, function(y) { 
    # Order by census date
    y_ord <- y[order(y[[census_date]]),]

    # Calculate stem growth
    y_ord[[w]][2] - y_ord[[w]][1]
  }))

  # Return
  return(out)
}

