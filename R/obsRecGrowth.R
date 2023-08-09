#' Estimate growth of individuals which recruited between two censues
#' 
#' `r descrip_table()` `r descrip_gro("calculate", "growth", "recruited")`
#'
#' @param x `r param_x()`
#' @param t0 `r param_t0()`
#' @param tT `r param_tT()`
#' @param w `r param_w()`
#' @param group `r param_group()`
#' @param census `r param_census()`
#' @param rec_method `r param_rec_method()`
#' @param diam optional, `r param_diam()`
#' @param min_size_class optional, `r param_min_size_class()`
#' @param min_diam_thresh optional, `r param_min_diam_thresh()`
#' @param growth_percentile optional, `r param_growth_percentile()`
#' @param w_min_diam optional, `r param_w_min_diam()`
#'
#' @return 
#' `r details_obs_sum()` growth from individuals which recruited.
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
#' obsRecGrowth(bicuar, "2019", "2021", w = "diam",
#'   group = "stem_id", census = "census_date",
#'   rec_method = "zero")
#' 
#' obsRecGrowth(bicuar, "2019", "2021", w = "diam",
#'   group = "stem_id", census = "census_date",
#'   rec_method = "thresh", w_min_diam = 5)
#' 
#' bicuar$agb_min <- runif(nrow(bicuar))
#' obsRecGrowth(bicuar, "2019", "2021", w = "agb",
#'   group = "stem_id", census = "census_date",
#'   rec_method = "thresh", w_min_diam = "agb_min")
#' 
#' obsRecGrowth(bicuar, "2019", "2021", w = "diam",
#'   group = "stem_id", census = "census_date",
#'   rec_method = "extrap", diam = "diam", min_size_class = c(5, 10),
#'   min_diam_thresh = 5, growth_percentile = 0.86)
#' 
#' @importFrom stats quantile
#' 
#' @export
#' 
obsRecGrowth <- function(x, t0, tT, w, group, census, 
  rec_method = "zero", diam = NULL, min_size_class = NULL, 
  min_diam_thresh = NULL, growth_percentile = NULL, w_min_diam = NULL) { 

  # Stop if methods not recognised
  if (any(!rec_method %in% c("zero", "thresh", "extrap"))) {
    stop("Invalid recruit growth estimation method. Methods must be either 'zero', 'thresh', 'extrap'")
  }

  # Stop if arguments required for method not provided, or malformed
  if (rec_method == "extrap") {
    if (is.null(diam) | is.null(min_size_class) | is.null(min_diam_thresh) | 
      is.null(growth_percentile)) {
      stop("Arguments required for rec_method == 'extrap' not provided")
    } 
    if (length(min_size_class) != 2 | !is.numeric(min_size_class)) {
      stop("min_size_class must be a numeric vector of length two")
    }
  }

  if (rec_method == "thresh") {
    if (is.null(w_min_diam)) {
      stop("Arguments required for rec_method == 'thresh' not provided")
    }
    # If w_min_diam is numeric, add as column
    if (is.numeric(w_min_diam)) {
      if (length(w_min_diam) > 1) {
        stop("Numeric w_min_diam must be of length one")
      }
      x$w_min_diam <- w_min_diam
      w_min_diam <- "w_min_diam"
    } else {
      if (!w_min_diam %in% names(x)) {
        stop("w_min_diam must be a numeric vector or a column name")
      }
    }
  }

  # Stop if any columns not recognised
  if (any(!c(w, group, census) %in% names(x))) {
    stop("Some columns not present in x")
  }

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

  # Find recruits
  ri <- obsID(x_fil, type = "rec", group = group, census = census,
      t0 = t0, tT = tT)

  # Create list to hold output
  method_list <- vector("list", length(rec_method))
  names(method_list) <- rec_method

  # If there are recruits:
  if (length(ri) > 0) {
    # Filter to recruits
    x_ri <- merge(x_fil, ri)

    # Filter recruits to final census
    x_ri_fil <- x_ri[x_ri[[census]] == tT,]

    # Method: Assume grew from 0 since previous census (zero)
    if ("zero" %in% rec_method) {
      # Calculate sum of growth
      BD <- sum(x_ri_fil[[w]], na.rm = TRUE)

      # Create list
      method_list[["zero"]] <- BD
    }

    # Method: Assume grew from min. diam. thresh. 
    if ("thresh" %in% rec_method) {
      # Find change
      BD <- sum(x_ri_fil[[w]] - x_ri_fil[[w_min_diam]], na.rm = TRUE)

      # Create list
      method_list[["thresh"]] <- BD
    }

    # Extrap. diam. to first census using percentile of small stem growth rate (extrap)
    if ("extrap" %in% rec_method) {
      # Find survivors
      si <- obsID(x_fil, type = "sur", group = group, census = census, 
        t0 = t0, tT = tT)

      # Filter to survivors
      x_si <- merge(x_fil, si)

      # Find stems in smallest diameter size class
      x_si_min <- x_si[x_si[[diam]] >= min_size_class[1] & 
        x_si[[diam]] < min_size_class[2],]

      # Find percentile growth rate of survivors in smallest diam size class 
      si_diam_growth_median <- stats::quantile(
        indGrowth(x_si_min, w = diam, group = group, census = census,
          t0 = t0, tT = tT), 
        growth_percentile)

      # For each stem, extrapolate back the stem diameter 
      x_ri_fil$diam_cen1 <- x_ri_fil$diam - si_diam_growth_median

      # Filter to extrapolated diameters > min. diam. thresh
      x_ri_gemin <- x_ri_fil[x_ri_fil$diam_cen1 > min_diam_thresh,]

      # Find change
      BD <- sum(x_ri_gemin[[w]], na.rm = TRUE)

      # Create list
      method_list[["extrap"]] <- BD
    }
  # If no recruits
  } else {
    # Fill with zero
    method_list[seq_along(method_list)] <- 0
  }

  # If only one method, don't return nested list
  if (length(rec_method) == 1) {
    method_list <- method_list[[1]]
  }

  # Return
  return(method_list)
}

