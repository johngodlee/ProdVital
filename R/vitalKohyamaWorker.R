#' Worker function to estimate vital rates according to 
#'    Kohyama et al. (2018), for single census interval in single plot
#'
#' @param x dataframe of SEOSAW stem data from single plot
#' @param w value from which to calculate vital rates, defaults to number of 
#'     stems if \code{NULL}
#' @param ind_id column name of individual IDs 
#' @param census_date column name of census dates 
#' @param plot_area plot area in hectares
#' @param census_date_1 column initial census of interval
#' @param census_date_2 column final census of interval
#' 
#' @return named list containing measures of vital rates
#' 
#' @export
#' 
vitalKohyamaWorker <- function(x, w = NULL, ind_id = "stem_id", 
  census_date = "census_date", plot_area, census_date_1, census_date_2) {

  # Recode plot area
  A <- plot_area 

  # Calculate census interval
  int <- as.numeric(census_date_2) - as.numeric(census_date_1)

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
  x_fil <- x[x[[census_date]] %in% c(census_date_1, census_date_2),]

  # If default measure of abundance
  if (is.null(w)) {
    x_fil$n <- 1
  } else {
    x_fil$n <- x_fil[[w]]
  }

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

  # Find initial abundance of living stems for plot
  N0 <- sum(x0$n, na.rm = TRUE)

  # Find final abundance of living stems for plot
  NT <- sum(xT$n, na.rm = TRUE)

  # Find stem IDs survivors
  si <- obsSur(x_fil, ind_id = ind_id, census_date = census_date, 
    census_date_1 = census_date_1, census_date_2 = census_date_2)

  # Find abundance of survivors
  NST <- sum(xT[xT[[ind_id]] %in% si, "n"], na.rm = TRUE)

  # Find stem IDs of recruits
  ri <- obsRec(x_fil, census_date = census_date, 
    census_date_1 = census_date_1, census_date_2 = census_date_2)

  # Find abundance of recruits
  nT <- sum(xT[xT[[ind_id]] %in% ri, "n"], na.rm = TRUE)

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

  # Calculate per area annual mortality rate
  Ma <- (N0 / A) * (1 - (NST / N0)^(1/int))

  # Calculate per area annual recruitment rate
  Ra <- Ma * (NT - NST) / (N0 - NST)

  # Calculate per area annual recruitment rate with first year deaths
  Ras <- Ma * (NT - NST) * (N0 / NST)^(1/int) / (N0 - NST)

  # Calculate per area annual recruitment rate with "Gf-estimate"
  #RGf <- G * (min_diam_thresh) * N * (min_diam_thresh) / A

  # Create output list
  out <- list(
    int = int,
    t0 = census_date_1,
    tT = census_date_2,
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

