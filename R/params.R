param_x <- function() {
  "dataframe of SEOSAW stem data"
}

param_w <- function() {
  "column name of value from which to calculate growth in `x`"
}

param_t0 <- function() {
  "initial census of interval"
}

param_tT <- function() {
  "final census of interval"
}

param_group <- function() {
  "column name(s) defining individual IDs in `x`"
}

param_census <- function() {
  "column name of census dates in `x`"
}

param_diam <- function() {
  "column name of stem diameters in `x`"
}

param_rec_method <- function() {
  "character string describing method used to estimate growth of recruits. Can be `zero`, `thresh`, or `extrap`."
}

param_growth_percentile <- function() {
  "percentile of growth rate used to estimate growth rate of recruits for backward extrapolation"
}

param_min_size_class <- function() {
  "vector of length two containing range of minimum diameter size class used to estimate median growth rate"
}

param_min_diam_thresh <- function() {
  "minimum diameter threshold measured"
}

param_w_min_diam <- function() {
  "column name or a vector of length 1 defining the estimated value of `w` in `x` at the minimum diameter threshold of the plot" 
}

param_size_class <- function() {
  "numeric vector containing diameter size class cut points. The largest number will be open ended, e.g. 50+."
}

param_plot_area <- function() {
  "plot area in hectares"
}

param_type <- function() {
  "vector defining which individuals to return, either: `rec` for recruits, `sur` for survivors, or `mor` for deaths"
}

param_type2 <- function() {
  "string defining whether to return only consecutive or all pairwise census intervals."
}

param_full <- function() {
  "logical, if `TRUE` return a dataframe, otherwise a named vector"
}

details_rec_method <- function() {
  "If `rec_method = 'thresh'`, `w_min_diam` must be defined. If `rec_method = 'extrap'`, `diam`, `min_size_class`, `min_diam_thresh` and `growth_percentile` must be defined."
}

details_growth_percentile <- function() {
  "Talbot et al. (2014) recommend a value of 0.86 for `growth_percentile` if `rec_method = 'extrap'`."
}

details_group <- function() {
  "`group` can be a vector of multiple character strings which, when combined, uniquely define each individual."
}

details_obs_sum <- function(un = FALSE) {
  if (un) {
    un_string <- "un"
  } else { 
    un_string <- ""
  }

  paste0("Named vector or dataframe (if `full = TRUE`) of individual ", 
    un_string, "observed")
}

descrip_table <- function() {
  "This function acts on a datafame (`x`) of individuals (`group`) measured at two time points (`t0`, `tT`, `census`)"
}

descrip_gro <- function(x = "calculate", y = "growth", z = "survived") {
  paste0("to ", x, " the ", y, " from individuals which ", z, ", given a measure of individual size (`w`)")
}

details_ncensus <- function() {
  "Only individuals with >1 census are returned."
}

