#' Impute missing and inconsistent mortality status
#'
#' @param x `r param_x()`
#' @param group `r param_group()`
#' @param census `r param_census()`
#' @param status column name of census dates in `x`
#' 
#' @return numeric vector of length \code{nrow(x)} with estimated individual 
#'     mortality status (1 = alive, 0 = dead)
#' 
#' @details 
#' Missing mortality status values are interpolated where possible. If the
#' missing value(s) is both preceded and followed by the same mortality status,
#' this value is used to fill the missing value(s). If the missing value(s) is
#' the first or final census(es) and either directly preceded or followed by two or more
#' of the same mortality status, this value is used to fill the missing
#' value(s). After these steps, all values occurring before the final alive
#' value are changed to alive. `r details_group()`
#' 
#' @examples
#' status_list <- list(
#'   c(NA_real_, NA_real_, 0, 0),
#'   c(NA_real_, NA_real_, 1, 0),
#'   c(NA_real_, NA_real_, 1, 1),
#'   c(NA_real_, NA_real_, 1),
#'   c(1, NA_real_, 1, 0),
#'   c(1, NA_real_, NA_real_, 0),
#'   c(1, NA_real_, NA_real_, 0),
#'   c(1, 0, 1, 0),
#'   c(1, NA_real_, 1, NA_real_))
#' 
#' lapply(status_list, function(x) { 
#'   d <- data.frame(
#'     stem_id = "a", 
#'     census_date = seq(2000, length.out = length(x)),
#'     status = x)
#' 
#'   mortFill(d, group = "stem_id",
#'     census = "census_date", status = "status")
#' })
#' 
#' @export
#' 
mortFill <- function(x, group, census, status) {

  # Create ID column
  x$row_id <- seq_len(nrow(x))

  # Order by census date
  x_ord <- x[order(x[[census]]),]

  # Generate unique individual IDs
  xint <- interaction(x[,group])

  # Split by group
  x_split <- split(x_ord, xint, drop = TRUE)

  # Impute dodgy and missing mortality values 
  x_fix <- do.call(rbind, lapply(x_split, function(y) { 
    # For each individual

    # If any status measurements are NA, but some are not
    if (any(is.na(y[[status]])) & any(!is.na(y[[status]]))) {

      # Compute run length encoding, excluding NAs
      y_rle_nona <- rle(c(na.omit(y[[status]])))

      # Compute run length encoding, with NAs
      y_rle <- rleNA(y[[status]])

      # Interpolate dead
      # If any dead found
      if (any(y_rle_nona$lengths[y_rle_nona$values == 0] > 0)) {
        # For each run of values
        for (i in seq_along(y_rle$values)[-1]) {
          # If status NA, preceded and followed by dead, not final value
          if (is.na(y_rle$values[i]) & 
            all(y_rle$values[c(i-1, i+1)] == 0, na.rm = TRUE) & 
            i != length(y_rle$values)) {
            # Fill NAs with alive 
            y_rle$values[i] <- 0
          }
        }
      }

      # Interpolate alive
      # If any alive found
      if (any(y_rle_nona$lengths[y_rle_nona$values == 1] > 0)) {
        # For each run of values
        for (i in seq_along(y_rle$values)[-1]) {
          # If status NA, preceded and followed by alive, not final value
          if (is.na(y_rle$values[i]) & 
            all(y_rle$values[c(i-1, i+1)] == 1, na.rm = TRUE) & 
            i != length(y_rle$values)) {
            # Fill NAs with alive 
            y_rle$values[i] <- 1
          }
        }
      }

      # If missing is first group
      if (is.na(y_rle$values[1])) {
        # If proceeded by 2 dead
        if (y_rle$values[2] == 0 & y_rle$lengths[2] >1) {
          # Fill NAs with dead
          y_rle$values[1] <- 0
        # If proceeded by 2 alive
        } else if (y_rle$values[2] == 1 & y_rle$lengths[2] >1) {
          # Fill NAs with alive
          y_rle$values[1] <- 1
        }
      }

      # If missing is last group
      if (is.na(y_rle$values[length(y_rle$values)])) {
        # If preceded by 2 dead
        if (y_rle$values[length(y_rle$values)-1] == 0 & 
          y_rle$lengths[length(y_rle$lengths)-1] >1) {
          # Fill NAs with dead
          y_rle$values[length(y_rle$values)] <- 0
        # If preceded by 2 alive
        } else if (y_rle$values[length(y_rle$values)-1] == 1 & 
          y_rle$lengths[length(y_rle$lengths)-1] >1) {
          # Fill NAs with alive
          y_rle$values[length(y_rle$values)] <- 1
        }
      }

      # Reverse run length encoding to get vector with interpolated values 
      y[[status]] <- inverse.rle(y_rle)
    }

    # Find alives
    alives <- which(y[[status]] == 1)

    # If any alives
    if (length(alives) > 0) {
      # Convert all within first and last to alive
      interm_change <- intersect(seq(alives[1], alives[length(alives)]),
        which(y[[status]] == 0))
      if (length(interm_change) > 0) {
        y[interm_change, status] <- 1
      }

      # Convert all before last alive to alive
      before_change <- 1:max(alives)
      if (length(before_change) > 0) {
        y[before_change, status] <- 1
      }
    }

    # Return dataframe with adjusted status and type of change
    data.frame(
      row_id = y$row_id,
      imput = y[[status]])
  }))

  # Order by ID column 
  out <- x_fix[order(x_fix$row_id), "imput"]

  # Return
  return(out)
}

