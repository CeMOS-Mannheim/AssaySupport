
#' Get normalization factors for a specific species
#'
#' @param peaksdf       Data.frame, see \code{AssaySupport::generate_peakdf}.
#' @param speciesdf     Data.frame, see \code{AssaySupport::generate_assigndf}.
#' @param targetSpecies Character, name of species for which the normalization factors should be extracted.
#' @param tol           Numeric, tolerance in Dalton or (if \code{tolppm = TRUE}) ppm
#' @param tolppm        Logical, use tolerance in ppm.
#' @param allowNoMatch  Logical, exclude spectra where \code{targetSpecies} was not found.
#'
#' @return
#' Numeric vector of normalization factors.
#' @export

getNormFactors <- function(peaksdf, speciesdf, targetSpecies, tol, tolppm = TRUE, allowNoMatch = TRUE) {
  
  resdf <- generate_resultdf(peaksdf, speciesdf = speciesdf, tolppm = TRUE, tol = tol)
  plot_Idx <- sort(unique(resdf$plotIdx))
  f_resdf <- resdf %>%
    filter(species == targetSpecies) %>% 
    arrange(plotIdx)
  
  
  
  if(!all(plot_Idx %in% (f_resdf %>% pull(plotIdx)))){
    if(!allowNoMatch) {
      stop("Could not find ", targetSpecies, " for all spectra! Consider adjusting tol.\n")
    }
    warning("Could not find ", targetSpecies, " in spectrum ", paste(which(!(plot_Idx %in% (f_resdf %>% pull(plotIdx)))), collapse = ", "), ".\n")
    specIdx <- which(plot_Idx %in% (f_resdf %>% pull(plotIdx)))
  } else {
    specIdx <- plot_Idx
  }
  return(list(norm_factor = pull(f_resdf, int),
              specIdx = specIdx))
}
