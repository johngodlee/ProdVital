#' Estimate vital rates according to Kohyama et al. (2018), across all
#'     census combinations for a single plot
#'
#' @param x dataframe of SEOSAW stem data from single plot
#' @param w value from which to calculate vital rates, defaults to number of 
#'     stems if \code{NULL}
#' @param ind_id column name of individual IDs 
#' @param census_date column name of census dates 
#' @param plot_area plot area in hectares
#' 
#' @return dataframe containing measures of vital rates for each pairwise 
#'     census interval:
#' \itemize{
#'   \item{int}{census interval length, in years}
#'   \item{t0}{year of initial census = census_date_1}
#'   \item{tT}{year of final census = census_date_2}
#'   \item{N0}{Abundance of living stems in initial census}
#'   \item{NT}{Abundance of living stems in final census}
#'   \item{NST}{Abundance of surviving stems in final census}
#'   \item{nT}{Abundance of recruits in final census}
#'   \item{g}{Intrinsic rate of natural increase}
#'   \item{m}{Instantaneous per capita mortality rate}
#'   \item{r}{Instantaneous per capita recruitment rate}
#'   \item{l}{Finite rate of increase}
#'   \item{ma}{Per capita annual mortality rate}
#'   \item{ra}{Initial-density-based per capita annual recruitment rate}
#'   \item{raf}{Final-density-based per capita annual recruitment rate}
#'   \item{ras}{Survivor-density-based per capita annual recruitment rate}
#'   \item{raz}{Zero-mortality per capita annual recruitment rate}
#'   \item{M}{Instantaneous per area mortality rate}
#'   \item{R}{Instantaneous per area recruitment rate}
#'   \item{Ma}{Per area annual mortality rate}
#'   \item{Ra}{Per area annual recruitment rate}
#'   \item{Ras}{Per area annual recruitment rate with first year deaths}
#' }
#' 
#' @export
#' 
vitalKohyama <- function(x, w = NULL, ind_id = "stem_id", 
  census_date = "census_date", plot_area) {

  # Convert potential tibble to dataframe
  x <- as.data.frame(x)

  # Find all census dates
  census_date_all <- sort(unique(x[[census_date]]))

  # Stop if only one census
  if (length(census_date_all) < 2) { 
    stop("Plot has only one census")
  }

  # Create all pairwise combinations of censuses
  comb_list_pair <- combn(census_date_all, 2, simplify = FALSE)

  # For each pairwise census interval combination:
  vital_pair_list <- lapply(comb_list_pair, function(i) { 
    # Estimate vital rates according to Kohyama 
    vitalKohyamaWorker(x, w = w, ind_id = ind_id, 
      census_date = census_date, plot_area = plot_area,
      census_date_1 = i[1], census_date_2 = i[2])
  })

  # Create dataframe of metrics
  out <- data.frame()[seq_along(vital_pair_list), ]
  for (i in names(vital_pair_list[[1]])) {
    out[[i]] <- unlist(lapply(vital_pair_list, "[[", i))
  }

  # Return
  return(out)
}

