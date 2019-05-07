#' Relaxed value matching with multiple results
#'
#' @param x   numeric vector or value, value(s) to match
#' @param vec numeric vector, to match \code{x} against
#' @param tol numeric, tolerance for matching
#' @param n   numeric, maximal number of matches. Set to "all" for all matches within \code{tol}
#'
#' @return
#' Integer vector, indeces of \code{vec} matching \code{x}.
#' @export

match_mutiple <- function(x, vec, tol = Inf, n = 1) {
  if(n == "all") {
    n <- length(vec)
  }
  tmp_dist <- dist <- abs(vec - x)
  res_vec <- vector("integer", length = n)
  for(i in 1:n) {
    idx <- which.min(dist)
    if(dist[idx] <= tol) {
      res_vec[i] <- idx
      dist[res_vec[i]] <- dist[res_vec[i]] + tol * 2
    }
  }
  return(res_vec)
}