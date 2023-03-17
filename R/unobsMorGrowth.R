#' Estimate growth of individuals which died between two censues
#'
#' @param x dataframe of stem measurements
#' @param w column name of value from which to calculate growth in \code{x}
#' @param ind_id column name of ind IDs in code{x}
#' @param diam column name of stem diameters in code{x}
#' @param census_date column name of census dates in code{x}
#' @param census_date_1 column initial census of interval
#' @param census_date_2 column final census of interval
#' @param size_class numeric vector containing diameter size class cut points.
#'     The largest number will be open ended, e.g. 50+
#'
#' @return total productivity from unobserved growth of stems which died
#' 
#' @details These are stems which were observed as alive in the previous 
#'     census, but then died sometime before the next census.
#' @export
#' 
unobsMorGrowth <- function(x, w = "diam", ind_id = "stem_id", diam = "diam", 
  census_date = "census_date", size_class = c(5,10,20,30,40,50),
  census_date_1, census_date_2) { 

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
    c(ind_id, w, census_date, diam)]

  # Find survivors between two censuses
  si <- obsSur(x_fil, ind_id = ind_id, census_date = census_date, 
    census_date_1 = census_date_1, census_date_2 = census_date_2)
  x_si <- x_fil[x_fil[[ind_id]] %in% si, 
    c(ind_id, w, census_date, diam)]

  # Define diameter size classes for each survivor in first census
  x_si_cen1 <- x_si[x_si[[census_date]] == census_date_1,]
  x_si_cen1$size_class <- as.character(cut(x_si_cen1[[diam]], breaks = size_class, 
    include.lowest = TRUE, right = FALSE, labels = size_class[-length(size_class)]))
  x_si_cen1$size_class[x_si_cen1[[diam]] > size_class[length(size_class)]] <- 
    paste0(size_class[length(size_class)], "+")
  x_si_cen1 <- x_si_cen1[,c(ind_id, "size_class")]

  # Apply size class to original dataframe
  x_si_class <- merge(x_si, x_si_cen1, 
    by.x = ind_id, by.y = ind_id, all.x = TRUE)

  # Define diameter size classes for each dead in first census
  x_di_cen1 <- x_di[x_di[[census_date]] == census_date_1,]
  x_di_cen1$size_class <- as.character(cut(x_di_cen1[[diam]], breaks = size_class, 
    include.lowest = TRUE, right = FALSE, labels = size_class[-length(size_class)]))
  x_di_cen1$size_class[x_di_cen1[[diam]] > size_class[length(size_class)]] <- 
    paste0(size_class[length(size_class)], "+")
  x_di_cen1 <- x_di_cen1[,c(ind_id, "size_class")]

  # Apply size class to original dataframe
  x_di_class <- merge(x_di, x_di_cen1, 
    by.x = ind_id, by.y = ind_id, all.x = TRUE)

  # Split surivors by size class
  x_si_split <- split(x_si_class, x_si_class$size_class)

  # Calculate median growth rate of survivors in all size classes
  si_class_median_growth <- lapply(x_si_split, function(y) {
    median(
      indGrowth(y, w = w, ind_id = ind_id, census_date = census_date, 
        census_date_1 = census_date_1, census_date_2 = census_date_2), 
      na.rm = TRUE)
  })

  # Create nice dataframe
  si_class_median_growth_df <- data.frame(
    size_class = names(si_class_median_growth),
    median_growth = unname(unlist(si_class_median_growth)))

  # Filter dead to first census
  x_di_class_cen1 <- x_di_class[x_di_class[[census_date]] == census_date_1,]

  # Join size class median growth rates to original dataframe 
  x_di_growth <- merge(x_di_class_cen1, si_class_median_growth_df, 
    by.x = "size_class", by.y = "size_class", all.x = TRUE)

  # Calculate diam growth assuming die half way (*0.5) through census interval
  x_di_growth$w_new <- x_di_growth[[w]] + (x_di_growth$median_growth * 0.5)

  # Create dataframe to estimate change
  x_cen1 <- x_di_growth[,c(ind_id, w)]
  x_cen2 <- x_di_growth[,c(ind_id, "w_new")]
  names(x_cen1) <- c(ind_id, w)
  names(x_cen2) <- c(ind_id, w)

  # Calculate woody productivity, 
  out <- sum(x_cen2[[w]] - x_cen1[[w]])

  # Return
  return(out)
}

