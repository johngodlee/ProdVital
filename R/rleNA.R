#' Relative length encoding, treating NA as a real value
#'
#' @param x vector
#'
#' @return relative length encoding of vector
#' 
#' @examples
#' x <- c(0,0,0,NA,NA,NA,1,1,1)
#' rle(x)
#' rleNA(x)
#' 
#' @export
#' 
rleNA <- function (x)  {
  stopifnot("'x' must be a vector of an atomic type" = is.atomic(x))

  n <- length(x)
  if (n == 0L) {
    return(structure(list(
      lengths = integer(), values = x)
    ), class = 'rle')
  }
  y <- (x[-1L] != x[-n])
  y[is.na(y)] <- FALSE
  y <- y | xor(is.na(x[-1L]), is.na(x[-n]))
  i <- c(which(y), n)

  structure(list(
    lengths = diff(c(0L, i)),
    values  = x[i]
  ), class = 'rle')
}


