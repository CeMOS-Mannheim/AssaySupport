#' Normalize spectra by factor
#'
#' @param spec     list, list of \code{MassSpectrum} objects.
#' @param factors  numeric vector, scaling factors. Needs to have same length as \code{spec}.
#'
#' @return
#' Scaled spectra
#' @export
normalizeByFactor <- function(spec, factors) {
  if(!length(spec) == length(factors)) {
    stop("Number of spectra and normalization factors not equal!\n")
  }
  
  spec_res <- spec
  
  for(i in 1:length(factors)) {
    spec_res[[i]]@intensity <- intensity(spec_res[[i]])/factors[i]
  }
  return(spec_res)
}

