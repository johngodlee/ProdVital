#' Estimate productivity according to Kohyama et al. (2019)
#'
#' `r descrip_table()` to estimate rates of productivity and loss.
#' 
#' @param x `r param_x()` from single plot
#' @param w `r param_w()`
#' @param t0 `r param_t0()`
#' @param tT `r param_tT()`
#' @param group `r param_group()`
#' @param census `r param_census()`
#' 
#' @return named list containing measures of productivity for the census interval: 
#'   * `int` - census interval length
#'   * `t0` - year of initial census 
#'   * `tT` - year of final census 
#'   * `N0` - number of living stems in initial census
#'   * `NT` - number of living stems in final census
#'   * `Ns0` - number of surviving stems in initial census
#'   * `Nr0` - number of recruiting stems 
#'   * `Nd0` - number of stems which died
#'   * `B0` - initial biomass
#'   * `BT` - final biomass
#'   * `Bs0` - biomass of survivors in initial census
#'   * `Br0` - biomass of recruits in final census
#'   * `Bd0` - biomass lost to deaths 
#'   * `dB` - observed net biomass change `BT - B0`
#'   * `dB_ann` - observed net annual biomass change `dB / int`
#'   * `W_max` - standardised maximum stem biomass for initial census (99th percentile of biomass)
#'   * `r_turn` - instantaneous recruitment rate (uses `turnoverEst()`)
#'   * `m_turn` - instantaneous mortality rate (uses `turnoverEst()`)
#'   * `p_turn` - instantaneous production rate (uses `turnoverEst()`)
#'   * `l_turn` - instantaneous loss rate 
#'   * `Bw` - period mean biomass `(BT-B0)/log(BT/B0)`
#'   * `Nw` - period mean abundance `(NT-N0)/log(NT/N0)`
#'   * `Bw_ann` - annual mean biomass `(BT-B0)/((BT/B0)^(1/int) - 1)/int`
#'   * `Nw_ann` - annual mean abundance `(NT-N0)/((NT/N0)^(1/int) - 1)/int`
#'   * `P` - instantaneous rate of production `(Bw * log(BT/Bs0))/int`
#'   * `L` - instantaneous rate of loss `(Bw * log(B0/Bs0))/int`
#'   * `P_simple` - simple rate of production `(B0-Bs0)/int`
#'   * `L_simple` - simple rate of loss `(BT-Bs0)/int`
#'   * `p` - intrinsic rate of production `(log(BT/Bs0)/int`
#'   * `l` - intrinsic rate of loss `(log(B0/Bs0)/int`
#'   * `r` - intrinsic rate of biomass change `(log(BT/B0)/int`
#'   * `p_ann` - specific rate of production `(BT/B0)^(1/int) * (1 - (Bs0/BT)^(1/int))`
#'   * `l_ann` - specific rate of loss `1-(Bs0/B0)^(1/int)`
#'   * `lambda` - specific rate of biomass change `(BT/B0)^(1/int)`
#'   * `P_ann` - annual rate of production `Bw_ann * p_ann`
#'   * `L_ann` - annual rate of loss `Bw_ann * l_ann`
#'   * `Psimp` - simple rate of production, alternative method `(NT*BT - Ns0*B0) / int`
#'   * `Psimp_clark` - simple rate of production (Clark et al. 2001) `Ns0*(BT-B0) / int + Nr0*(BT-Bmin) / int`
#' 
#' @details
#' `r details_group()`
#'
#' @examples
#' data(bicuar_clean)
#' bicuar_p1 <- bicuar_clean[bicuar_clean$plot_id == "ABG_5",]
#' 
#' prodKohyama(bicuar_p1, "2019", "2021", w = "diam",
#'   group = "stem_id", census = "census_date")
#' 
#' @importFrom stats quantile
#' 
#' @export
#' 
prodKohyama <- function(x, t0, tT, w, group, census) {

  # Calculate census interval
  int <- as.numeric(tT) - as.numeric(t0)

  # Stop is interval is negative
  if (int < 0) { 
    stop("census dates must be in chronological order")
  }

  # Subset to censuses of interest
  x_fil <- x[x[[census]] %in% c(t0, tT),]
  x0 <- x_fil[x_fil[[census]] == t0,]
  xT <- x_fil[x_fil[[census]] == tT,]

  # If no stems found in a census, throw warning
  if (nrow(x0) == 0) {
    no1 <- TRUE
    warning("Census 1 contains no stems") 
  } 
  if (nrow(xT) == 0) {
    no2 <- TRUE
    warning("Census 2 contains no stems")
  }

  # Define minimum biomass possible on plot
  Bmin <- 0

  # Find initial number of stems for plot
  N0 <- nrow(x0)

  # Find final number of stems for plot
  NT <- nrow(xT)

  # Find initial biomass for plot
  B0 <- sum(x0[[w]], na.rm = TRUE)

  # Find final biomass for plot
  BT <- sum(xT[[w]], na.rm = TRUE)

  # Find change in biomass for plot
  dB <- BT - B0

  # Find annual change in biomass for plot
  dB_ann <- dB / int

  # Find survivors
  si <- obsID(x_fil, type = "sur", group = group, census = census, 
        t0 = t0, tT = tT)
  Ns0 <- nrow(si)
  
  # Find deaths
  di <- obsID(x_fil, type = "mor", group = group, census = census,
      t0 = t0, tT = tT)
  Nd0 <- nrow(di)

  # Find recruits
  ri <- obsID(x_fil, type = "rec", group = group, census = census,
      t0 = t0, tT = tT)
  Nr0 <- nrow(ri)

  # Should be equal in a perfect dataset
  n0 <- Ns0 + Nd0  # Initial N equals survivors plus deaths
  # n0 == N0
  nT <- Ns0 + Nr0  # Final N equals survivors plus recruits
  # nT == NT  # NT is sometimes larger than nT, 
  # implying some stems appear as alive in final census, 
  # but not found or are dead/resprouting in initial census.

  # Find initial biomass for survivors
  x_si <- merge(x_fil, si)
  Bs0 <- sum(x_si[x_si[[census]] == t0, w], na.rm = TRUE)

  # Find final biomass for recruits 
  x_ri <- merge(x_fil, ri)
  Br0 <- sum(x_ri[x_ri[[census]] == t0, w], na.rm = TRUE)

  # Find biomass lost to deaths 
  x_di <- merge(x_fil, di)
  Bd0 <- sum(x_di[x_di[[census]] == t0, w], na.rm = TRUE)

  # Find standardised maximum stem biomass for initial census
  W_max <- as.numeric(stats::quantile(x_fil[
    !interaction(x_fil[,group]) %in% interaction(ri) & 
      x_fil[[census]] == t0, 
    w], probs = 0.99, na.rm = TRUE))

  # Eq10
  # Period mean biomass 
  Bw <- ifelse(BT != B0, (BT-B0)/log(BT/B0), B0)
  # Period mean abundance 
  Nw <- ifelse(NT != N0, (NT-N0)/log(NT/N0), N0)

  # Eq14
  # Annual mean biomass of census period
  Bw_ann <- ifelse(BT != B0, (BT-B0)/((BT/B0)^(1/int) - 1)/int, B0)
  # Annual mean abundance of census period
  Nw_ann <- ifelse(NT != N0, (NT-N0)/((NT/N0)^(1/int) - 1)/int, N0)

  # Instantaneous rates
  P = (Bw * log(BT/Bs0))/int  # Production Eq1
  L = (Bw * log(B0/Bs0))/int  # Loss Eq2

  # Simple rates 
  P_simple  = (BT - Bs0)/int  # Production Eq5
  L_simple  = (B0 - Bs0)/int  # Loss Eq6

  # Intrinsic rates of biomass change
  p <- log(BT / Bs0) / int  # Production Eq8
  l <- log(B0 / Bs0) / int  # Loss Eq9
  r <- log(BT / B0) / int  # Net change Eq7-ish

  # Specific rates
  p_ann <- (BT/B0)^(1/int) * (1 - (Bs0/BT)^(1/int))  # Production Eq12
  l_ann <- 1 - (Bs0/B0)^(1/int)  # Loss Eq13
  lambda <- (BT/B0)^(1/int)  # Net change (Lambda)

  # Annual rates
  P_ann = Bw_ann * p_ann  # Production Eq3
  L_ann = Bw_ann * l_ann  # Loss Eq4

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

  Psimp <- (NT * BT - Ns0 * B0) / int
  Psimp_clark  <- Ns0 * (BT - B0) / int + Nr0 * (BT - Bmin) / int

  Tw <- B0 / P_ann

  # Create output list
  out <- list(
    int = int,
    t0 = t0,
    tT = tT,
    N0 = N0, NT = NT, Ns0 = Ns0, Nr0 = Nr0, Nd0 = Nd0,
    B0 = B0, BT = BT, Bs0 = Bs0, Br0 = Br0, Bd0 = Bd0,
    dB = dB, dB_ann = dB_ann,
    W_max = W_max,
    r_turn = r_turn, m_turn = m_turn, p_turn = p_turn, l_turn = l_turn,
    Bw = Bw, Nw = Nw, 
    Bw_ann = Bw_ann, Nw_ann = Nw_ann,
    P = P, L = L,
    P_ann = P_ann, L_ann = L_ann, 
    P_simple = P_simple, L_simple = L_simple, 
    p = p, l = l, r = r,
    p_ann = p_ann, l_ann = l_ann, lambda = lambda,
    Psimp = Psimp, Psimp_clark = Psimp_clark,
    Tw = Tw
)

  # Return
  return(out)
}

