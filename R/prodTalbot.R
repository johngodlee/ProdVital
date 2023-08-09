#' Estimate productivity according to Talbot et al. (2014)
#' 
#' `r descrip_table()` to estimate rates of producitivity.
#'
#' @param x `r param_x()` from single plot
#' @param t0 `r param_t0()`
#' @param tT `r param_tT()`
#' @param w `r param_w()`
#' @param diam `r param_diam()`
#' @param group `r param_group()`
#' @param census `r param_census()`
#' @param size_class `r param_size_class()` Passed to `unobsMorGrowth()`
#' @param min_size_class `r param_min_size_class()`
#' @param w_min_diam `r param_w_min_diam()`
#' @param rec_method `r param_rec_method()`
#' @param min_diam_thresh optional, `r param_min_diam_thresh()`
#' @param growth_percentile optional, `r param_growth_percentile()`
#'
#' @return named list containing measures of productivity:
#'   * `t0` - census date of initial census
#'   * `tT` - census date of final census
#'   * `int` - `tT - t0`
#'   * `sur_obs` - Observed growth of survivors
#'   * `rec_obs` - Observed growth contribution of recruits observed in final census
#'   * `mor_obs` - Observed loss from stems which died in census interval
#'   * `mor_unobs` - Estimated unobserved growth of stems which died in census interval
#'   * `rec_unobs` - Estimated unobserved growth of stems which recruited and died in census interval
#'   * `AGWP_obs` - Observed productivity: `sur_obs + rec_obs`
#'   * `AGWP_est` - Estimated productivity using CIC2, stem by stem approach: `AGWP_obs + mor_unobs + rec_unobs`
#'   * `AGWP_obs_ann` - `AGWP_obs / int`
#'   * `AGWP_est_ann` - `AGWP_est / int`
#' 
#' @details 
#' `r details_group()`
#' 
#' `r details_rec_method()`
#' 
#' `r details_growth_percentile()`
#'
#' @examples
#' data(bicuar)
#' 
#' bicuar$agb_min <- runif(nrow(bicuar))
#' prodTalbot(bicuar, "2019", "2021", w = "agb", diam = "diam",
#'   group = "stem_id", census = "census_date", w_min_diam = "agb_min", 
#'   rec_method = "zero")
#' 
#' prodTalbot(bicuar, "2019", "2021", w = "agb", diam = "diam",
#'   group = "stem_id", census = "census_date", w_min_diam = "agb_min",
#'   size_class = seq(5, 20, 5), min_size_class = c(5, 10),
#'   rec_method = "zero")
#' 
#' prodTalbot(bicuar, "2019", "2021", w = "diam", diam = "diam",
#'   group = "stem_id", census = "census_date", w_min_diam = 5, 
#'   rec_method = "thresh")
#' 
#' prodTalbot(bicuar, "2019", "2021", w = "agb", diam = "diam",
#'   group = "stem_id", census = "census_date", w_min_diam = "agb_min",
#'   rec_method = "thresh")
#' 
#' prodTalbot(bicuar, "2019", "2021", w = "agb", diam = "diam", 
#'   group = "stem_id", census = "census_date", w_min_diam = "agb_min", 
#'   rec_method = "extrap", min_diam_thresh = 5, growth_percentile = 0.86)
#'   
#' @export
#' 
prodTalbot <- function(x, t0, tT, w, diam, group, census, 
  size_class = c(5,10,20,30,40,50), min_size_class = c(5, 10), 
  w_min_diam, rec_method = "zero", 
  min_diam_thresh = NULL, growth_percentile = NULL) {

  # Calculate census interval
  int <- as.numeric(tT) - as.numeric(t0)

  # Stop is interval is negative
  if (int < 0) { 
    stop("census dates must be in chronological order")
  }

  # Subset to censuses of interest
  x_fil <- x[x[[census]] %in% c(t0, tT),]

  # If no stems found in a census, throw warning
  if (!t0 %in% x_fil[[census]]) {
    warning("Census 1 contains no stems") 
  }
  if (!tT %in% x_fil[[census]]) {
    warning("Census 2 contains no stems") 
  }

  # Calculate growth of survivors
  obs_sur_growth <- obsSurGrowth(x_fil, w = w, group = group, 
    census = census,
    t0 = t0, tT = tT)

  # Calculate growth of observed recruits
  obs_rec_growth <- obsRecGrowth(x_fil, w = w, group = group, 
    diam = diam, census = census,
    growth_percentile = growth_percentile,
    min_size_class = min_size_class,
    rec_method = rec_method,
    min_diam_thresh = min_diam_thresh,
    t0 = t0, tT = tT,
    w_min_diam = w_min_diam)

  # Calculate unobserved growth of stems which died
  unobs_mor_growth <- unobsMorGrowth(x_fil, w = w, group = group, 
    diam = diam, census = census,
    size_class = size_class,
    t0 = t0, tT = tT)

  # Calculate unobserved growth of unobserved recruits
  unobs_rec_growth <- unobsRecGrowth(x_fil, t0 = t0, tT = tT, 
    w = w, group = group, census = census, diam = diam, 
    min_size_class = min_size_class, w_min_diam = w_min_diam)

  # Calculate observed loss to mortality
  obs_mor_loss <- obsMorLoss(x_fil, w = w, group = group, 
    census = census, 
    t0 = t0, tT = tT)

  # Total change = growth of survivors + growth of recruits
  AGWP_obs <- obs_sur_growth + obs_rec_growth

  # Estimate total woody productivity = observed change + estimates from dead stems
  AGWP_est <- obs_sur_growth + obs_rec_growth + 
    unobs_mor_growth + unobs_rec_growth

  # Annual rates
  AGWP_est_ann <- AGWP_est / int
  AGWP_obs_ann <- AGWP_obs / int

  # Create list of values
  out <- list(t0 = t0, tT = tT, int = int, 
    sur_obs = obs_sur_growth, rec_obs = obs_rec_growth,
    mor_obs = obs_mor_loss, 
    mor_unobs = unobs_mor_growth, rec_unobs = unobs_rec_growth,
    AGWP_obs = AGWP_obs, AGWP_est = AGWP_est, 
    AGWP_obs_ann = AGWP_obs_ann, AGWP_est_ann = AGWP_est_ann)

  # Return
  return(out)
}

