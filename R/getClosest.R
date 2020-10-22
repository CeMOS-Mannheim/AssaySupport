#' get closest matches
#'
#' @param x       Numeric vector
#' @param y       Numeric vector to match x to.
#' @param n       Numeric, maximum number of matches
#' @param tol     Numeric, tolerance for matching.
#' @param unique  Logical, only retrun unique matches (if n>1).
#'
#' @export
getClosest <- function (x, y, n = NA, tol = NA, unique = FALSE) 
{
  diff1 <- diff2 <- abs(x - y)
  if (is.na(n) && is.na(tol) && !unique) {
    return(which(diff1 %in% min(diff1)))
  }
  if (!is.na(tol)) {
    diff1[which(!diff1 < tol)] <- Inf
    if (all(diff1 == Inf)) {
      return(NA)
    }
  }
  res <- which(diff1 %in% sort(diff2, partial = n)[1:n])
  res <- res[order(diff1[res])]
  if (unique) {
    res <- res[!duplicated(x[res])]
  }
  return(res)
}
