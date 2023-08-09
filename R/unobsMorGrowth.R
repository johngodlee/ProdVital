#' Estimate growth of individuals which died between two censues
#' 
#' `r descrip_table()` `r descrip_gro("estimate", "growth", "died")`
#'
#' @param x `r param_x()`
#' @param t0 `r param_t0()`
#' @param tT `r param_tT()`
#' @param w `r param_w()`
#' @param group `r param_group()`
#' @param census `r param_census()`
#' @param diam `r param_diam()`
#' @param size_class `r param_size_class()`
#'
#' @return 
#' `r details_obs_sum(un = TRUE)` growth from individuals which died.
#' 
#' @details 
#' `r details_group()`
#'
#' @examples
#' data(bicuar)
#' 
#' unobsMorGrowth(bicuar, "2019", "2021", w = "diam", 
#'   group = "stem_id", census = "census_date", diam = "diam")
#' 
#' @importFrom stats median
#' 
#' @export
#' 
unobsMorGrowth <- function(x, t0, tT, w, group, census,
  diam, size_class = c(5,10,20,30,40,50)) { 

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
  x_di <- merge(x_fil, di)

  # Find survivors between two censuses
  si <- obsID(x_fil, type = "sur", group = group, census = census, 
    t0 = t0, tT = tT)
  x_si <- merge(x_fil, si)

  # Define diameter size classes for each survivor in first census
  x_si_cen1 <- x_si[x_si[[census]] == t0,]
  x_si_cen1$size_class <- as.character(cut(x_si_cen1[[diam]], breaks = size_class, 
    include.lowest = TRUE, right = FALSE, labels = size_class[-length(size_class)]))
  x_si_cen1$size_class[x_si_cen1[[diam]] > size_class[length(size_class)]] <- 
    paste0(size_class[length(size_class)], "+")
  x_si_cen1 <- x_si_cen1[,c(group, "size_class")]

  # Apply size class to original dataframe
  x_si_class <- merge(x_si, x_si_cen1, 
    by.x = group, by.y = group, all.x = TRUE)

  # Define diameter size classes for each dead in first census
  x_di_cen1 <- x_di[x_di[[census]] == t0,]
  x_di_cen1$size_class <- as.character(cut(x_di_cen1[[diam]], breaks = size_class, 
    include.lowest = TRUE, right = FALSE, labels = size_class[-length(size_class)]))
  x_di_cen1$size_class[x_di_cen1[[diam]] > size_class[length(size_class)]] <- 
    paste0(size_class[length(size_class)], "+")
  x_di_cen1 <- x_di_cen1[,c(group, "size_class")]

  # Apply size class to original dataframe
  x_di_class <- merge(x_di, x_di_cen1, 
    by.x = group, by.y = group, all.x = TRUE)

  # Split surivors by size class
  x_si_split <- split(x_si_class, x_si_class$size_class)

  # Calculate median growth rate of survivors in all size classes
  si_class_median_growth <- lapply(x_si_split, function(y) {
    stats::median(
      indGrowth(y, w = w, group = group, census = census, 
        t0 = t0, tT = tT), 
      na.rm = TRUE)
  })

  # Create nice dataframe
  si_class_median_growth_df <- data.frame(
    size_class = names(si_class_median_growth),
    median_growth = unname(unlist(si_class_median_growth)))

  # Create unrepresented size classes
  unrep_size_class <- unique(x_di_class$size_class)[
    !unique(x_di_class$size_class) %in% si_class_median_growth_df$size_class]

  if (length(unrep_size_class) > 0) {
    if (nrow(si_class_median_growth_df) > 0) {
      median_growth_all <- mean(si_class_median_growth_df$median_growth, 
        na.rm = TRUE)
    } else {
      median_growth_all <- 0
    }

    unrep_size_class_df <- data.frame(
      size_class = unrep_size_class,
      median_growth = median_growth_all)

    si_class_median_growth_df <- rbind(
      si_class_median_growth_df, unrep_size_class_df)
  }

  # Filter dead to first census
  x_di_class_cen1 <- x_di_class[x_di_class[[census]] == t0,]

  # Join size class median growth rates to original dataframe 
  x_di_growth <- merge(x_di_class_cen1, si_class_median_growth_df, 
    by.x = "size_class", by.y = "size_class", all.x = TRUE)

  # Calculate diam growth assuming die half way (*0.5) through census interval
  x_di_growth$w_new <- x_di_growth[[w]] + (x_di_growth$median_growth * 0.5)

  # Create dataframe to estimate change
  x_cen1 <- x_di_growth[,c(group, w)]
  x_cen2 <- x_di_growth[,c(group, "w_new")]
  names(x_cen1) <- c(group, w)
  names(x_cen2) <- c(group, w)

  # Calculate woody productivity, 
  out <- sum(x_cen2[[w]] - x_cen1[[w]], na.rm = TRUE)

  # Return
  return(out)
}

