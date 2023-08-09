#' Find individuals which survived, died, or recruited 
#' 
#' `r descrip_table()` to return the identities of individuals which survived
#' (measured at both time points), recruited (measured only at the final time
#' point), or died (measured only at the initial time point).
#'
#' @param x `r param_x()`
#' @param t0 `r param_t0()`
#' @param tT `r param_tT()`
#' @param type `r param_type()`
#' @param group `r param_group()`
#' @param census `r param_census()`
#' 
#' @details 
#' `r details_group()`
#' 
#' `r details_t0_tT()`
#'
#' @return 
#' Dataframe or list of dataframes of IDs for individuals which recruited,
#' survived, or died, between censuses. 
#' 
#' @examples
#' data(bicuar)
#' 
#' obsID(bicuar, "2019", "2021", 
#'   group = "stem_id", census = "census_date")
#' 
#' obsID(bicuar, "2019", "2021", type = "mor", 
#'   group = "stem_id", census = "census_date")
#' 
#' @export
#' 
obsID <- function(x, t0, tT, type = c("rec", "sur", "mor"), group, census) {

  # Check validity of type column
  if (!all(type %in% c("rec", "sur", "mor"))) {
    stop("'type' must contain only values of 'rec', 'sur', and 'mor'")
  }

  # Convert potential tibble to dataframe
  x <- as.data.frame(x)

  # Subset to initial census
  x0 <- x[x[[census]] == t0, group, drop = FALSE]

  # Subset to final census
  xt <- x[x[[census]] == tT, group, drop = FALSE]

  # Generate unique individual IDs
  x0int <- interaction(x0[,group])
  xtint <- interaction(xt[,group])

  # Find stem IDs in x0 but not in xt
  outl <- list()
  for (i in seq_along(type)) {
    if (type[i] == "rec") {
      outl[[i]] <- xt[!xtint %in% x0int, , drop = FALSE]
    } else if (type[i] == "mor") {
      outl[[i]] <- x0[!x0int %in% xtint, , drop = FALSE]
    } else if (type[i] == "sur") {
      outl[[i]] <- x0[x0int %in% xtint, , drop = FALSE]
    }
  }
  names(outl) <- type

  # Simplify if only one item in type
  if (length(outl) == 1) {
    outl <- outl[[1]]
  }

  # Return
  return(outl)
}


