#' Worker function to estimate woody productivity according to 
#'     Talbot et al. (2014), for single census interval in single plot
#'
#' @param x dataframe of SEOSAW stem data from single plot
#' @param w column name of value from which to calculate growth in \code{x}
#' @param diam column name of stem diameters in code{x}
#' @param ind_id column name of stem IDs in code{x}
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
#' @param census_date_1 column initial census of interval
#' @param census_date_2 column final census of interval
#' @param w_min the minimum value of \code{w} expected to be encountered in 
#'     the plot given the minimum diameter threshold. Default assumes minimum 
#'     diameter threshold
#'
#' @return named list containing measures of productivity 
#' 
#' @export
#' 
prodTalbotWorker <- function(x, w = "agb", diam = "diam", 
  ind_id = "stem_id", census_date = "census_date", 
  min_diam_thresh = 5, growth_percentile = 0.86, min_size_class = c(5,10),
  rec_method = "min_diam_thresh", size_class = c(5,10,20,30,40,50),
  census_date_1, census_date_2, w_min = min_size_class[1]) {

  # Calculate census interval
  int <- as.numeric(census_date_2) - as.numeric(census_date_1)

  # Stop is interval is negative
  if (int < 0) { 
    stop("census dates must be in chronological order")
  }

  # Subset to censuses of interest
  x_fil <- x[x[[census_date]] %in% c(census_date_1, census_date_2),]

  # If no stems found in a census, throw warning
  if (!census_date_1 %in% x_fil[[census_date]]) {
    warning("Census 1 contains no stems") 
  }
  if (!census_date_2 %in% x_fil[[census_date]]) {
    warning("Census 2 contains no stems") 
  }

  # Calculate growth of survivors
  obs_sur_growth <- obsSurGrowth(x_fil, w = w, ind_id = ind_id, 
    census_date = census_date,
    census_date_1 = census_date_1, census_date_2 = census_date_2)

  # Calculate growth of observed recruits
  obs_rec_growth <- obsRecGrowth(x_fil, w = w, ind_id = ind_id, 
    diam = diam, census_date = census_date,
    growth_percentile = growth_percentile,
    min_size_class = min_size_class,
    rec_method = rec_method,
    min_diam_thresh = min_diam_thresh,
    census_date_1 = census_date_1, census_date_2 = census_date_2,
    w_min = w_min)

  # Calculate unobserved growth of stems which died
  unobs_mor_growth <- unobsMorGrowth(x_fil, w = w, ind_id = ind_id, 
    diam = diam, census_date = census_date,
    size_class = size_class,
    census_date_1 = census_date_1, census_date_2 = census_date_2)

  # Calculate unobserved growth of unobserved recruits
  unobs_rec_growth <- unobsRecGrowth(x_fil, w = w, ind_id = ind_id,
    diam = diam, census_date = census_date,
    min_size_class = min_size_class,
    census_date_1 = census_date_1, census_date_2 = census_date_2, 
    w_min = w_min)

  # Total change = growth of survivors + growth of recruits
  AGWP_obs <- obs_sur_growth + obs_rec_growth

  # Estimate total woody productivity = observed change + 
  AGWP_est <- obs_sur_growth + obs_rec_growth + 
    unobs_mor_growth + unobs_rec_growth

  # Annual rates
  AGWP_est_ann <- AGWP_est / int
  AGWP_obs_ann <- AGWP_est / int

  # Create list of values
  out <- list(AGWP_obs = AGWP_obs, AGWP_est = AGWP_est, 
    AGWP_obs_ann = AGWP_obs_ann, AGWP_est_ann = AGWP_est_ann, 
    sur_obs = obs_sur_growth, rec_obs = obs_rec_growth,
    mor_unobs = unobs_mor_growth, rec_unobs = unobs_rec_growth,
    t0 = census_date_1, tT = census_date_2, int = int)

  # Return
  return(out)
}

