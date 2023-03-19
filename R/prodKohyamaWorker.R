#' Worker function to estimate woody productivity according to 
#'    Kohyama et al. (2019), for single census interval in single plot
#'
#' @param x dataframe of SEOSAW stem data from single plot
#' @param ind_id column name of individual IDs 
#' @param agb column name of AGB values
#' @param census_date column name of census dates 
#' @param census_date_1 column initial census of interval
#' @param census_date_2 column final census of interval
#' 
#' @return named list containing measures of productivity
#' 
#' @export
#' 
prodKohyamaWorker <- function(x, ind_id = "stem_id", 
  agb = "agb", census_date = "census_date", census_date_1, census_date_2) {

  # Calculate census interval
  int <- as.numeric(census_date_2) - as.numeric(census_date_1)

  # Stop is interval is negative
  if (int < 0) { 
    stop("census dates must be in chronological order")
  }

  # Subset to censuses of interest
  x_fil <- x[x[[census_date]] %in% c(census_date_1, census_date_2),]
  x0 <- x_fil[x_fil[[census_date]] == census_date_1,]
  xT <- x_fil[x_fil[[census_date]] == census_date_2,]

  # If no stems found in a census, throw warning
  if (nrow(x0) == 0) {
    no1 <- TRUE
    warning("Census 1 contains no stems") 
  } 
  if (nrow(xT) == 0) {
    no2 <- TRUE
    warning("Census 2 contains no stems")
  }

  # Find initial number of stems for plot
  N0 <- nrow(x0)

  # Find final number of stems for plot
  NT <- nrow(xT)

  # Find initial biomass for plot
  B0 <- sum(x0[[agb]], na.rm = TRUE)

  # Find final biomass for plot
  BT <- sum(xT[[agb]], na.rm = TRUE)

  # Find change in biomass for plot
  dB <- BT - B0

  # Find annual change in biomass for plot
  dB_ann <- dB / int

  # Find survivors
  si <- obsSur(x_fil, ind_id = ind_id, census_date = census_date, 
        census_date_1 = census_date_1, census_date_2 = census_date_2)
  Ns0 <- length(si)
  
  # Find deaths
  di <- obsMor(x_fil, ind_id = ind_id, census_date = census_date,
      census_date_1 = census_date_1, census_date_2 = census_date_2)
  Nd0 <- length(di)

  # Find recruits
  ri <- obsRec(x_fil, ind_id = ind_id, census_date = census_date,
      census_date_1 = census_date_1, census_date_2 = census_date_2)
  Nr0 <- length(ri)

  # Should be equal in a perfect dataset
  n0 <- Ns0 + Nd0  # Initial N equals survivors plus deaths
  # n0 == N0
  nT <- Ns0 + Nr0  # Final N equals survivors plus recruits
  # nT == NT  # NT is sometimes larger than nT, 
    # implying some stems appear as alive in final census, 
    # but not found or are dead/resprouting in initial census.

  # Find initial biomass for survivors
  Bs0 <- sum(x_fil[x_fil[[ind_id]] %in% si & 
      x_fil[[census_date]] == census_date_1, 
    agb], na.rm = TRUE)

  # Find final biomass for recruits 
  Br0 <- sum(x_fil[x_fil[[ind_id]] %in% ri & 
      x_fil[[census_date]] == census_date_2, 
    agb], na.rm = TRUE)

  # Find biomass lost to deaths 
  Bd0 <- sum(x_fil[x_fil[[ind_id]] %in% di & 
      x_fil[[census_date]] == census_date_1, 
    agb], na.rm = TRUE)

  # Find standardised maximum stem biomass for initial census
  W_max <- as.numeric(quantile(x_fil[
      !x_fil[[ind_id]] %in% ri & 
        x_fil[[census_date]] == census_date_1,
      agb], probs = 0.99, na.rm = TRUE))

  # Find minimum biomass found on plot 
  Bmin <- 0

  # Period mean abundance - Kohyama et al. (2018) Eq10
  Nw <- ifelse(NT != N0, (NT-N0)/log(NT/N0), N0)

  # Period mean biomass - Kohyama et al. (2018) Eq10
  Bwk <- ifelse(BT != B0, (BT-B0)/log(BT/B0), B0)
  
  # Kohyama et al. (2018) Eq14
  Nw_ann <- ifelse(NT != N0, (NT-N0)/((NT/N0)^(1/int) - 1)/int, N0)
  Bw_ann <- ifelse(BT != B0, (BT-B0)/((BT/B0)^(1/int) - 1)/int, B0)

  # Calculate production, loss, and mean biomass 
  L_simple  = sum(B0 - Bs0)/int  # Simple loss
  P_simple  = sum(BT - Bs0)/int  # Simple production
  Bw_simple = sum(B0 + BT)/2  # Simple mean biomass

  L_ann  = sum(Bw_ann * (1 - (Bs0/B0)^(1/int)))  # Annual loss
  P_ann  = sum(Bw_ann * (BT/B0)^(1/int) * (1 - (Bs0/BT)^(1/int)))  # Annual production

  L = sum(Bwk * log(B0/Bs0))/int  # Instantaneous loss
  P = sum(Bwk * log(BT/Bs0))/int  # Instantaneous production
  Bw = sum(Bwk)  # Period mean biomass

  # Calculate turnover rates
  r_turn <- try(turnoverEst(NT, Ns0, int), silent = TRUE)
  if (inherits(r_turn,"try-error")) {
    r_turn <- NA_real_
  }
  m_turn <- try(turnoverEst(N0, Ns0, int), silent = TRUE)
  if (inherits(m_turn,"try-error")) {
    m_turn <- NA_real_
  }
  p_turn <- try(turnoverEst(BT, Bs0, int), silent = TRUE)
  if (inherits(p_turn,"try-error")) {
    p_turn <- NA_real_
  }
  l_turn <- try(turnoverEst(B0, Bs0, int), silent = TRUE)
  if (inherits(l_turn,"try-error")) {
    l_turn <- NA_real_
  }

  Psimp <- (sum((NT * BT - Ns0 * B0) / int)) 
  Psimp_clark  <- (sum(Ns0 * (BT - B0) / int + Nr0 * (BT - Bmin) / int))

  # Create output list
  out <- list(
    int = int,
    t0 = census_date_1,
    tT = census_date_2,
    N0 = N0, NT = NT, Ns0 = Ns0, Nr0 = Nr0, Nd0 = Nd0,
    B0 = B0, BT = BT, Bs0 = Bs0, Br0 = Br0, Bd0 = Bd0,
    dB = dB, dB_ann = dB_ann,
    W_max = W_max,
    r_turn = r_turn, m_turn = m_turn, p_turn = p_turn, l_turn = l_turn,
    Nw = Nw, Nw_ann = Nw_ann,
    Bwk = Bwk, Bw_ann = Bw_ann,
    P_simple = P_simple, L_simple = L_simple, Bw_simple = Bw_simple,
    P_ann = P_ann, L_ann = L_ann,
    P = P, L = L, Bw = Bw,
    Psimp = Psimp, Psimp_clark = Psimp_clark)

  # Return
  return(out)
}

