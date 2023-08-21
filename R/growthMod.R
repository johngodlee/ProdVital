#' Estimate average growth rate of individuals
#'
#' @param x `r param_x()`
#' @param w `r param_w()`
#' @param group `r param_group()`
#' @param census `r param_census()`
#' @param full `r param_full()`
#'
#' @return 
#' `r details_obs_sum()` average growth rates across all census
#' intervals.
#'
#' @details
#' This function fits a linear model to estimate average growth rate, and
#' returns the slope of that model, which is identical to the average rate
#' of growth, in units of `w` per `census`. Only individuals with >1 census
#' are returned.
#' 
#' @examples
#' data(bicuar_clean)
#' 
#' growthMod(bicuar_clean, w = "diam",
#'   group = c("plot_id", "stem_id"), census = "census_date")
#' 
#' growthMod(bicuar_clean, w = "diam",
#'   group = c("plot_id", "stem_id"), census = "census_date", 
#'   full = TRUE)
#' 
#' @export
#' 
growthMod <- function(x, w, group, census, full = FALSE) {

  # Convert potential tibble to dataframe
  x <- as.data.frame(x)
 
  # Find individuals with more than one census
  ind_id <- apply(x[,group], 1, paste, collapse = ":")
  ind_count <- table(ind_id)
  ind_multi <- names(ind_count)[ind_count > 1]
  x_fil <- x[which(ind_id %in% ind_multi),]

  # Split dataframe by indvidual
  x_split <- split(x_fil, x_fil[,group], drop = TRUE)

  # Fit linear models and extract slope coefficient
  coefs <- unlist(lapply(x_split, function(i) {
    coef(.lm.fit(cbind(1, i[[census]]), i[[w]]))[2]
  }))

  # Add names or return dataframe
  if (full) {
    out <- unique(x_fil[,group])
    out$growth_rate <- coefs
  } else {
    out <- coefs
    names(out) <- unique(interaction(x_fil[,group]))
  }

  # Return
  return(out)
}

