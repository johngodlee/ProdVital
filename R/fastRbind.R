#' Fast rbind for dataframes
#'
#' @param x list of dataframe objects
#'
#' @return dataframe
#' 
#' @examples
#' x <- list(
#'   data.frame(x = seq(1, 5), y = seq(5, 9)),
#'   data.frame(x = seq(10, 50), y = seq(50, 90)))
#' fastRbind(x)
#' 
#' @export
#' 
fastRbind <- function(x) { 
  list2DF(lapply(setNames( seq_along(x[[1]]), names(x[[1]])), 
      function(i) {
        unlist(lapply(x, `[[`, i), FALSE, FALSE)
      }))
}

