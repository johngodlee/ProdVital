#' Calculate stem growth increment across all censuses in a plot
#' 
#' @param x dataframe of SEOSAW stem data 
#' @param w column name of value from which to calculate growth in \code{x}
#' @param ind_id column name of individual IDs 
#' @param census_date column name of census dates 
#' 
#' @return vector length \code{nrow(x)} with stem growth increments 
#'     for each census.
#' 
#' @export
#' 
indGrowthInc <- function(x, w = "diam", ind_id = "stem_id",
  census_date = "census_date") {

  # Convert to clean dataframe
  x <- as.data.frame(x[,c(w, ind_id, census_date)])

  # Create row ID for sorting
  x$row_id <- seq_len(nrow(x))

  # Extract all censuses
  census_year_all <- sort(unique(x[[census_date]]))

  # Create stepwise combinations
  start_date <- c(NA, census_year_all)
  end_date <- c(census_year_all, NA)

  # Calculate growth increment
  inc_df <- fastRbind(lapply(seq(2, length(start_date)-1), function(y) {
    growth_vec <- indGrowth(x = x, w = w, ind_id = ind_id, 
      census_date = census_date, 
      census_date_1 = start_date[y], census_date_2 = end_date[y])
    
    data.frame(
      ind_id = names(growth_vec),
      growth_inc = growth_vec,
      census_date = end_date[y])
  }))

  # Merge dataframes
  inc_merge <- merge(x, inc_df, 
    by.x = c(ind_id, census_date), by.y = c("ind_id", "census_date"), 
    all.x = TRUE)

  # Order by row ID
  out <- inc_merge[order(inc_merge$row_id), "growth_inc"]

  # Return
  return(out)
}

