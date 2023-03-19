#' Estimate woody productivity according to Talbot et al. (2014), across all
#'     census combinations for a single plot
#' 
#' @param x dataframe of SEOSAW stem data from single plot
#' @param w column name of value from which to calculate growth in \code{x}
#' @param ind_id column name of individual IDs in code{x}
#' @param diam column name of stem diameters in code{x}
#' @param census_date column name of census dates in code{x}
#' @param min_diam_thresh minimum diameter threshold of plot
#' @param growth_percentile percentile of growth rate used to estimate growth 
#'     rate of recruits for backward extrapolation
#' @param min_size_class vector of length two containing range of minimum 
#'     diameter size class used to estimate median growth rate of recruits
#' @param rec_method character string describing method used to estimate 
#'     growth of recruits
#' @param size_class numeric vector containing diameter size class cut points.
#'     The largest number will be open ended, e.g. 50+
#' @param w_min the minimum value of \code{w} expected to be encountered in 
#'     the plot given the minimum diameter threshold. Default assumes minimum 
#'     diameter threshold
#'
#' @return dataframe containing measures of prodcutivity for each pairwise census interval:
#' \itemize{
#'   \item{AGWP_obs}{Observed productivity: sur_obs + rec_obs}
#'   \item{AGWP_est}{Estimated productivity using CIC2, stem by stem approach: AGWP_obs + mor_unobs + rec_unobs}
#'   \item{AGWP_obs_ann}{AGWP_obs / int}
#'   \item{AGWP_est_ann}{AGWP_est / int}
#'   \item{sur_obs}{Observed growth of survivors}
#'   \item{rec_obs}{Observed growth contribution of recruits observed in final census}
#'   \item{mor_unobs}{Estimated unobserved growth of stems which died in census interval}
#'   \item{rec_unobs}{Estimated unobserved growth of stems which recruited and died in census interval}
#'   \item{mor_obs}{Observed loss from stems which died in census interval}
#'   \item{t0}{census date of initial census}
#'   \item{tT}{census date of final census}
#'   \item{int}{tT - t0}
#'   \item{cic}{correction factor calculated as slope of linear model of mean AGWP_obs and int}
#'   \item{AGWP_est_ann_cic}{Estimated productivity using CIC1, plot correction factor approach: AGWP_obs * (cic - int)}
#' }
#' 
#' @details Only returns pairwise estimates for census intervals 
#'     of less than 10 years.
#' @export
#' 
prodTalbot <- function(x, w = "agb", ind_id = "stem_id", diam = "diam",
  census_date = "census_date", 
  min_diam_thresh = 5, growth_percentile = 0.86, min_size_class = c(5,10),
  rec_method = "min_diam_thresh", size_class = c(5,10,20,30,40,50), w_min = min_size_class[1]) { 

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

  # Discard censuses >10 years apart
  comb_list_pair_fil <- comb_list_pair[unlist(lapply(comb_list_pair, diff)) <=10]

  # If no censuses <= 10 years apart, return nothing
  if (length(comb_list_pair_fil) > 0) {
    # For each pairwise census interval combination
    agwp_pair_list <- lapply(comb_list_pair_fil, function(i) { 
      prodTalbotWorker(x, w = w, diam = diam, ind_id = ind_id,
        census_date = census_date, 
        min_diam_thresh = min_diam_thresh, 
        growth_percentile = growth_percentile, min_size_class = min_size_class, 
        rec_method = rec_method, size_class = size_class, 
        census_date_1 = i[1], census_date_2 = i[2], w_min = w_min)
    })

    # Define window function
    window_func <- function(x, y) {
      lapply(seq_len(length(x) - y + 1), function(i) {
        x[i:(i + y - 1)]
      })
    }

    # Define window lengths
    window_len <- min(length(census_date_all), 4)
    window_seq <- seq(2, window_len)

    # Find all census interval combinations at 2,3,4 combinations
    comb_list <- unlist(lapply(window_seq, window_func, x = census_date_all), 
      recursive = FALSE)

    # Find all single interval combinations
    comb_list_two <- comb_list[unlist(lapply(comb_list, length)) == 2]

    # For each two census interval combination
    agwp_int_list <- lapply(comb_list_two, function(i) { 
      prodTalbotWorker(x, w = w, diam = diam, ind_id = ind_id,
        census_date = census_date, min_diam_thresh = min_diam_thresh, 
        growth_percentile = growth_percentile, min_size_class = min_size_class, 
        rec_method = rec_method, size_class = size_class, 
        census_date_1 = i[1], census_date_2 = i[2], w_min = w_min)
    })

    # Extract annual AGWP from each interval combination
    AGWP_obs_ann_vec <- unlist(lapply(agwp_int_list, "[[", "AGWP_obs_ann"))
    t0_vec <- unlist(lapply(agwp_int_list, "[[", "t0"))
    tT_vec <- unlist(lapply(agwp_int_list, "[[", "tT"))

    # For each combination calculate mean AGWP_ann from stepwise combinations
    comb_two_join <- unlist(lapply(comb_list_two, paste, collapse = "-"))
    comb_list_join <- lapply(comb_list, function(y) { 
      date_join <- paste(c(NA, y), c(y, NA), sep = "-") 
      date_join[-c(1, length(date_join))]
    })

    # For each census interval combination extract AGWP
    mean_agwp_df <- fastRbind(lapply(comb_list_join, function(y) {
      mean_AGWP_ann <- mean(AGWP_obs_ann_vec[which(comb_two_join %in% y)])
      mean_int <- mean(tT_vec[which(comb_two_join %in% y)] - 
        t0_vec[which(comb_two_join %in% y)])
      data.frame(mean_AGWP_ann, mean_int)
    }))

    # Run linear model of mean interval vs. mean annual AGWP
    cic_lm <- lm(mean_AGWP_ann ~ mean_int, data = mean_agwp_df)

    # Extract slope from model
    cic <- unname(cic_lm$coefficients[2])
    cic <- ifelse(is.na(cic), 0, cic)

    # Extract change for all pairwise census intervals
    out <- data.frame()[seq_along(agwp_pair_list), ]
    for (i in names(agwp_pair_list[[1]])) {
      out[[i]] <- unlist(lapply(agwp_pair_list, "[[", i))
    }
    out$cic <- cic
    out$AGWP_est_ann_cic <- out$AGWP_obs_ann - (out$cic * out$int)
    
    # Return dataframe
    return(out)
  }
}

