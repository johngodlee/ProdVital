#' Calculate growth increment for consecutive or pairwise intervals
#'
#' @param x `r param_x()`
#' @param w `r param_w()`
#' @param group `r param_group()`
#' @param census `r param_census()`
#' @param type `r param_type2()`
#'
#' @return 
#' Dataframe of individual growth increment values across all
#' consecutive or pairwise census intervals.
#' 
#' Columns:
#' 
#' * `t0` - date of the start of the census interval
#' * `tT` - date of the end of the census interval
#' * `w0` - The value of `w` at the start of the census interval
#' * `wT` - The value of `w` at the end of the census interval
#' * `g` - Growth increment (`wT - w0`)
#' 
#' @details
#' This function calculates growth increments between census
#' intervals. Column `t0` is the initial census in the census
#' interval, `tT` is the final census, `int` is the length of the
#' census, and `g` is the growth increment.
#' 
#' `type` must be one of "consecutive" or "pairwise".
#' 
#' `r details_ncensus()`
#' 
#' @examples
#' data(bicuar_clean)
#' 
#' growthAll(bicuar_clean, w = "diam",
#'   group = c("plot_id", "stem_id"), census = "census_date",
#'   type = "consecutive")
#' 
#' growthAll(bicuar_clean, w = "diam",
#'   group = c("plot_id", "stem_id"), census = "census_date",
#'   type = "pairwise")
#' 
#' @export
#' 
growthAll <- function(x, w, group, census, type = "consecutive") {

  # Convert potential tibble to dataframe
  x <- as.data.frame(x)

  # Check, is `type` specified correctly?
  if (!type %in% c("consecutive", "pairwise")) {
    stop("`type` must be one of 'consecutive' or 'pairwise'")
  }
 
  # Find individuals with more than one census
  ind_id <- apply(x[,group, drop = FALSE], 1, paste, collapse = ":")
  ind_count <- table(ind_id)
  ind_multi <- names(ind_count)[ind_count > 1]
  x_fil <- x[which(ind_id %in% ind_multi),]

  # Split by group
  x_split <- split(x_fil, x_fil[,group], drop = TRUE)

  # Calculate pairwise or consecutive growth increments
  out <- do.call(rbind, lapply(x_split, function(y) {
    # Find all census dates
    census_all <- sort(unique(y[[census]]))

    if (type == "pairwise") {
      # Create all pairwise combinations of censuses
      comb_list <- combn(census_all, 2, simplify = FALSE)
    } else if (type == "consecutive") {
      # Extract first and last censuses
      t0 <- c(NA_real_, census_all)
      tT <- c(census_all, NA_real_)

      # Create all consecutive combinations of censuses
      comb_list <- lapply(seq_along(t0), function(i) { 
        if (all(!is.na(c(t0[i], tT[i])))) {
          c(t0[i], tT[i]) 
        }
      })
      
      # Remove NULL entries
      comb_list <- comb_list[lengths(comb_list) > 0]
    }

    # For each in comb_list
    do.call(rbind, lapply(comb_list, function(i) { 
      # Find growth increment
      gi <- y[1, group, drop = FALSE]
      gi$t0 <- y[y[[census]] == i[1], census]
      gi$tT <- y[y[[census]] == i[2], census]
      gi$int <- gi$tT - gi$t0
      gi$w0 <- y[y[[census]] == i[1], w]
      gi$wT <- y[y[[census]] == i[2], w]
      gi$g <- gi$wT - gi$w0

      gi
    }))
  }))

  # Return
  return(out)
}


