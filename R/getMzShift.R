#' Get normalization factors for a specific species
#'
#' @param peaksdf       Data.frame, see \code{AssaySupport::generate_peakdf}.
#' @param speciesdf     Data.frame, see \code{AssaySupport::generate_assigndf}.
#' @param targetSpecies Character or Numeric, name of species for which the mz shift should be extracted.
#'                      Or, if numeric, use mz to match species.
#' @param tol           Numeric, tolerance in Dalton or (if \code{tolppm = TRUE}) ppm
#' @param tolppm        Logical, use tolerance in ppm.
#' @param allowNoMatch  Logical, exclude spectra where \code{targetSpecies} was not found.
#'
#' @return
#' List with two entries
#' mzshift: Numeric vector of mz shift per spectrum.
#' specIdx: Numeric vector of spectra indicies for which the targetSpecies was found.
#' 
#' @importFrom dplyr mutate arrange filter
#' 
#' @export
#' 

getMzShift <- function(peaksdf, tol, targetSpecies, speciesdf = NA, tolppm = TRUE, allowNoMatch = TRUE) {
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
  resdf <- generate_resultdf(peaksdf, speciesdf = speciesdf, tolppm = TRUE, tol = tol)
  f_resdf <- resdf %>%
    filter(species == targetSpecies) %>%
    group_by(plotIdx) %>%
    filter(abs(mz.diff) == min(abs(mz.diff))) %>%
    arrange(plotIdx)
  } else {
    if(tolppm) {
      resdf <- peaksdf %>%
        mutate(match = ifelse(mz > targetSpecies - mz*(tol/1e6) & mz < targetSpecies + mz*(tol/1e6), TRUE, FALSE)) 
      f_resdf <- resdf %>%
        filter(match) %>%
        mutate(mz.diff = round(targetSpecies - mz, 4)) %>%
        group_by(plotIdx) %>%
        filter(abs(mz.diff) == min(abs(mz.diff))) %>%
        arrange(plotIdx)
    } else {
  resdf <- peaksdf %>%
    mutate(match = mz > targetSpecies - tol & mz < targetSpecies + tol) 
  f_resdf <- resdf %>%
    filter(match) %>%
    mutate(mz.diff = round(targetSpecies - mz, 4)) %>%
    group_by(plotIdx) %>%
    filter(abs(mz.diff) == min(abs(mz.diff))) %>%
    arrange(plotIdx)
    }
  }
  
  if(!all(plot_Idx %in% (f_resdf %>% pull(plotIdx)))){
    if(!allowNoMatch) {
      stop("Could not find ", targetSpecies, " for all spectra! Consider adjusting tol.\n")
    }
    warning("Could not find ", targetSpecies, " in spectrum ", paste(which(!(plot_Idx %in% (f_resdf %>% pull(plotIdx)))), collapse = ", "), ".\n")
    specIdx <- sort(which(plot_Idx %in% (f_resdf %>% pull(plotIdx))))
  } else {
    specIdx <- plot_Idx
  }
  if(length(specIdx) < 1) {
    stop("Could not find targetSpecies in any spectrum! Consider adjusting tol.\n")
  }
  
  return(list(mzshift = pull(f_resdf, mz.diff),
              specIdx = specIdx))
}
