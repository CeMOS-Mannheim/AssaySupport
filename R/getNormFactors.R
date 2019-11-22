
#' Get normalization factors for a specific species
#'
#' @param peaksdf       Data.frame, see \code{AssaySupport::generate_peakdf}.
#' @param speciesdf     Data.frame, see \code{AssaySupport::generate_assigndf}.
#' @param targetSpecies Character or Numeric, name of species for which the normalization factors should be extracted.
#'                      Or, if numeric, use mz to match species.
#' @param tol           Numeric, tolerance in Dalton or (if \code{tolppm = TRUE}) ppm
#' @param tolppm        Logical, use tolerance in ppm.
#' @param allowNoMatch  Logical, exclude spectra where \code{targetSpecies} was not found.
#'
#' @return
#' List with two entries
#' norm_factor: Numeric vector of normalization factors per spectrum.
#' specIdx: Numeric vector of spectra indicies for which the targetSpecies was found
#' 
#' @examples 
#' plot(test_spectra_proc[[1]], xlim = c(3000,4000))
#' peakdf <- generate_peakdf(test_spectra_proc)
#' norm <- getNormFactors(peakdf, targetSpecies = 3316.4, tol = 100)
#' norm_spec <- normalizeByFactor(test_spectra_proc, norm$norm_factor)
#' plot(norm_spec[[1]], xlim = c(3000,4000))
#' 
#' @importFrom dplyr mutate arrange filter
#' 
#' @export

getNormFactors <- function(peaksdf, targetSpecies, tol, speciesdf = NA, tolppm = TRUE, allowNoMatch = TRUE) {
  
  if(is.character(targetSpecies)) {
    if(!is.data.frame(speciesdf)) {
      stop("If targetSpecies is specfied by a character string then speciesdf is needed.\n")
    }
    if(is.numeric(targetSpecies)) {
      if(is.data.frame(speciesdf)) {
        warning("speciesdf will not be used if targetSpecies is provided as numeric (mz value).\n")
      }
    }
  }
  
  plot_Idx <- sort(unique(peaksdf$plotIdx))
  
  
  if(is.character(targetSpecies)) {
    resdf <- generate_resultdf(peaksdf, speciesdf = speciesdf %>% filter(Species == targetSpecies), tolppm = TRUE, tol = tol)
    f_resdf <- resdf %>%
      filter(species == targetSpecies) %>%
      arrange(plotIdx)
  } else {
    if(tolppm) {
      resdf <- peaksdf %>%
        mutate(match = ifelse(mz > targetSpecies - mz*(tol/1e6) & mz < targetSpecies + mz*(tol/1e6), TRUE, FALSE)) 
      f_resdf <- resdf %>%
        filter(match) %>%
        arrange(plotIdx)
    } else {
      resdf <- peaksdf %>%
        mutate(match = mz > targetSpecies - tol & mz < targetSpecies + tol) 
      f_resdf <- resdf %>%
        filter(match) %>%
        arrange(plotIdx)
    }
  }
  
  
  
  if(!all(plot_Idx %in% (f_resdf %>% pull(plotIdx)))){
    if(!allowNoMatch) {
      stop("Could not find ", targetSpecies, " for all spectra! Consider adjusting tol.\n")
    }
    warning("Could not find ", targetSpecies, " in spectrum ", paste(which(!(plot_Idx %in% (f_resdf %>% pull(plotIdx)))), collapse = ", "), ".\n")
    specIdx <- which(plot_Idx %in% (f_resdf %>% pull(plotIdx)))
  } else {
    specIdx <- plot_Idx
  }
  if(length(specIdx) < 1) {
    stop("Could not find targetSpecies in any spectrum! Consider adjusting tol.\n")
  }
  return(list(norm_factor = pull(f_resdf, int),
              specIdx = specIdx))
}
