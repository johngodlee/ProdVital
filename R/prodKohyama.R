#' Estimate woody productivity according to Kohyama et al. (2019), across all
#'     census combinations for a single plot
#'   
#' @param x dataframe of SEOSAW stem data from single plot
#' @param ind_id column name of individual IDs 
#' @param agb column name of AGB values
#' @param census_date column name of census dates 
#' @param plot_area plot area in hectares
#'
#' @return dataframe containing measures of productivity for each pairwise 
#'     census interval:
#' \itemize{
#'   \item{int}{census interval length, in years}
#'   \item{t0}{year of initial census = census_date_1}
#'   \item{tT}{year of final census = census_date_2}
#'   \item{N0}{number of living stems in initial census}
#'   \item{NT}{number of living stems in final census}
#'   \item{Ns0}{number of surviving stems in initial census}
#'   \item{Nr0}{number of recruiting stems }
#'   \item{Nd0}{number of stems which died}
#'   \item{B0}{initial biomass}
#'   \item{BT}{final biomass}
#'   \item{Bs0}{biomass of survivors in initial census}
#'   \item{Br0}{biomass of recruits in final census}
#'   \item{Bd0}{biomass lost to deaths }
#'   \item{W_max}{standardised maximum stem biomass for initial census (99th percentile of biomass)}
#'   \item{r_turn}{instantaneous recruitment rate}
#'   \item{m_turn}{instantaneous mortality rate }
#'   \item{p_turn}{instantaneous production rate}
#'   \item{l_turn}{instantaneous loss rate }
#'   \item{Nw}{period mean abundance (Kohyama et al. 2019 - Eq10)}
#'   \item{Narea}{period mean abundance per hectare}
#'   \item{Nw_ann}{annual mean abundance (Kohyama et al. 2019 - Eq14) }
#'   \item{Bw}{period mean biomass (Kohyama et al. 2019 - Eq10)}
#'   \item{Barea}{period mean biomass per hectare}
#'   \item{Bw_ann}{annual mean biomass (Kohyama et al. 2019 - Eq14)}
#'   \item{P_simple}{simple rate of production (B0-Bs0)/int}
#'   \item{L_simple}{simple rate of loss (BT-Bs0)/int}
#'   \item{Bw_simple}{simple mean biomass (B0+BT)/2}
#'   \item{P_ann}{annual rate of production (Kohyama et al. 2019 Eq12)}
#'   \item{L_ann}{annual rate of loss (Kohyama et al. 2019 Eq13)}
#'   \item{P}{instantaneous rate of production}
#'   \item{L}{instantaneous rate of loss}
#'   \item{Bw}{alternative measure of period mean biomass, to check consistency in calculation with Bwk}
#'   \item{Pabs}{absolute productivity (p_turn * Barea)}
#'   \item{Psimp}{simple rate of production, alternative method}
#'   \item{Psimp_clark}{simple rate of production (Clark et al. 2001)}
#' }
#' 
#' @details Only returns pairwise estimates for census intervals
#'     of less than 10 years.
#' 
#' @export
#' 
prodKohyama <- function(x, ind_id = "stem_id", 
  agb = "agb", census_date = "census_date", plot_area) {

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

  # Discard censuses >10 years apart
  comb_list_pair_fil <- comb_list_pair[unlist(lapply(comb_list_pair, diff)) <=10]

  # If no censuses <= 10 years apart, return nothing
  if (length(comb_list_pair_fil) > 0) {
    # For each pairwise census interval combination:
    agwp_pair_list <- lapply(comb_list_pair_fil, function(i) { 
      # Estimate productivity according to Kohyama 
      prodKohyamaWorker(x, ind_id = ind_id, agb = agb, 
        census_date = census_date, plot_area, 
        census_date_1 = i[1], census_date_2 = i[2])
    })

    # Create dataframe of metrics
    out <- data.frame()[seq_along(agwp_pair_list), ]
    for (i in names(agwp_pair_list[[1]])) {
      out[[i]] <- unlist(lapply(agwp_pair_list, "[[", i))
    }

    # Return
    return(out)
  }
}

