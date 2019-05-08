#' Relaxed value matching with multiple results
#'
#' @param x       numeric, value to match
#' @param vec     numeric vector, to match \code{x} against
#' @param tol     numeric, tolerance for matching. 
#'                Multiply with 1e-6 and set \code{tol_ppm} to use relative tolerance.
#' @param n       numeric, maximal number of matches. Set to "all" for all matches within \code{tol}
#' @param tol_ppm logical, set to \code{TRUE} if \code{tol} should be used as ppm.
#'
#' @return
#' Integer vector, indeces of \code{vec} matching \code{x}.
#' @export

match_mutiple <- function(x, vec, tol = Inf, n = 1, tol_ppm = FALSE) {
  
  if(length(x) > 1) {
    stop("x has to be a single numeric value.\n")
  }
  
  if(n == "all") {
    n <- length(vec)
  } else if(!is.numeric(n)){
    stop("n has to be numeric or 'all'.\n")
  }
  
  dist <- abs(vec - x)
  res_vec <- vector("integer", length = n)
  res_vec[] <- NA
  for(i in 1:n) {
    idx <- which.min(dist)
    if(tol_ppm) {
      new_tol <- tol * x
    } else {
      new_tol = tol}
    if(dist[idx] <= new_tol) {
      res_vec[i] <- idx
      dist[res_vec[i]] <- dist[res_vec[i]] + tol * 1e6
    }
  }
  return(res_vec[which(!is.na(res_vec))])
}



