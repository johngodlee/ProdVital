#' Estimate vital rates according to Kohyama et al. (2018)
#' 
#' `r descrip_table()` to estimate rates of recruitment, mortality, and net
#' population change.
#'
#' @param x `r param_x()` from single plot
#' @param t0 `r param_t0()`
#' @param tT `r param_tT()`
#' @param w `r param_w()` defaults to number of stems if \code{NULL}
#' @param group `r param_group()`
#' @param census `r param_census()`
#' @param plot_area `r param_plot_area()`
#' 
#' @return named list containing measures of vital rates for census interval:
#'   * `int` - census interval length, in years
#'   * `t0` - year of initial census = t0
#'   * `tT` - year of final census = tT
#'   * `N0` - Abundance of living stems in initial census
#'   * `NT` - Abundance of living stems in final census
#'   * `NST` - Abundance of surviving stems in final census
#'   * `nT` - Abundance of recruits in final census
#'   * `g` - Intrinsic rate of natural increase
#'   * `m` - Instantaneous per capita mortality rate
#'   * `r` - Instantaneous per capita recruitment rate
#'   * `l` - Finite rate of increase
#'   * `ma` - Per capita annual mortality rate
#'   * `ra` - Initial-density-based per capita annual recruitment rate
#'   * `raf` - Final-density-based per capita annual recruitment rate
#'   * `ras` - Survivor-density-based per capita annual recruitment rate
#'   * `raz` - Zero-mortality per capita annual recruitment rate
#'   * `M` - Instantaneous per area mortality rate
#'   * `R` - Instantaneous per area recruitment rate
#'   * `Ma` - Per area annual mortality rate
#'   * `Ra` - Per area annual recruitment rate
#'   * `Ras` - Per area annual recruitment rate with first year deaths
#'
#' @details
#' `r details_group()`
#' 
#' @examples
#' data(bicuar)
#' 
#' vitalKohyama(bicuar, "2019", "2021", group = "stem_id", 
#'   census = "census_date", plot_area = 1)
#' 
#' vitalKohyama(bicuar, "2019", "2021", group = "stem_id", w = "agb",
#'   census = "census_date", plot_area = 1)
#' 
#' @export
#' 
vitalKohyama <- function(x, t0, tT, w = NULL, group, census, plot_area) {

  # Recode plot area
  A <- plot_area 

  # Calculate census interval
  int <- as.numeric(tT) - as.numeric(t0)

  # Stop is interval is negative
  if (int < 0) { 
    stop("census dates must be in chronological order")
  }

  # If alternative measure of abundance check that column is valid
  if (!is.null(w)) {
    if (is.null(x[[w]])) {
      stop("'w' is not a valid column in 'x'")
    }
  }

  # Subset to censuses of interest
  x_fil <- x[x[[census]] %in% c(t0, tT),]

  # If default measure of abundance
  if (is.null(w)) {
    x_fil$n <- 1
  } else {
    x_fil$n <- x_fil[[w]]
  }

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

  # Find initial abundance of living stems for plot
  N0 <- sum(x0$n, na.rm = TRUE)

  # Find final abundance of living stems for plot
  NT <- sum(xT$n, na.rm = TRUE)

  # Find stem IDs survivors
  si <- obsID(x_fil, type = "sur", group = group, census = census, 
    t0 = t0, tT = tT)
  xT_si <- merge(xT, si)

  # Find abundance of survivors
  NST <- sum(xT_si$n, na.rm = TRUE)

  # Find stem IDs of recruits
  ri <- obsID(x_fil, type = "rec", group = group, census = census, 
    t0 = t0, tT = tT)
  xT_ri <- merge(xT, ri)

  # Find abundance of recruits
  nT <- sum(xT_ri$n, na.rm = TRUE)

  # Calculate intrinsic rate of natural increase
  g <- log(NT / NST) / int

  # Calculate instantaneous per capita mortality rate
  m <- log(N0 / NST) / int

  # Calculate instantaneous per capita recruitment rate
  r <- log(NT / NST) / int

  # Calculate finite rate of increase
  l <- (NT / N0)^(1/int)

  # Calculate per capita annual mortality rate
  ma <- 1 - (NST / N0)^(1/int)

  # Calculate initial-density-based per capita annual recruitment rate
  ra <- ( (NT / NST)^(1/int) - 1) * ((NST / N0)^(1-int))

  # Calculate final-density-based per capita annual recruitment rate
  raf <- 1 - (NST / NT) ^ (1/int)

  # Calculate survivor-density-based per capita annual recruitment rate
  ras <- (NT / NST)^(1/int) - 1

  # Calculate zero-mortality per capita annual recruitment rate
  raz <- (1 + NT/N0 - NST/N0)^(1/int) - 1

  # Calculate instantaneous per area mortality rate
  M <- (N0 / A) * log(N0 / NST) / int

  # Calculate instantaneous per area recruitment rate
  R <- M * (NT - NST) / (N0 - NST)
  if (is.nan(R)) {
    R <- (NT / A) * log(NT / NST) / int
  }

  # Calculate per area annual mortality rate
  Ma <- (N0 / A) * (1 - (NST / N0)^(1/int))

  # Calculate per area annual recruitment rate
  Ra <- Ma * (NT - NST) / (N0 - NST)
  if (is.nan(Ra)) {
    Ra <- (NT / A) * (1 - (NST / NT)^(1/int))
  }

  # Calculate per area annual recruitment rate with first year deaths
  Ras <- Ma * (NT - NST) * (N0 / NST)^(1/int) / (N0 - NST)
  if (is.nan(Ras)) {
    Ras <- Ra
  }

  # Calculate per area annual recruitment rate with "Gf-estimate"
  #RGf <- G * (min_diam_thresh) * N * (min_diam_thresh) / A

  # Create output list
  out <- list(
    int = int,
    t0 = t0,
    tT = tT,
		N0 = N0,
		NT = NT,
		NST = NST,
		nT = nT,
		g = g,
		m = m,
		r = r,
		l = l,
		ma = ma,
		ra = ra,
		raf = raf,
		ras = ras,
		raz = raz,
		M = M,
		R = R,
		Ma = Ma,
		Ra = Ra,
		Ras = Ras)

  # Return
  return(out)
}

