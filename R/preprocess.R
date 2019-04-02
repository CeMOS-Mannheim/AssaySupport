#' Preprocess spectra
#'
#' @param spec                  List of \code{MALDIquant::MassSpectrum}
#' @param smooth_meth           Character, smoothing method see \code{MALDIquant::smoothIntensity()}
#' @param smooth_halfWindowSize Numeric, halfWindowSize for smoothing method see \code{MALDIquant::smoothIntensity()}
#' @param norm_meth             Character, normalization method see \code{MALDIquant::calibrateIntensity()}
#' @param filter_spectra        Character vector, regex of patterns to exclude spectra (like calibration spectra in same folder)
#' @param baseline_meth         Character, baseline removal method see \code{MALDIquant::removeBaseline}
#' @param avg_method            Character, aggregation method used to generate average spectra. See \code{MALDIquant::averageMassSpectra}.
#' @param align                 Chacter vector, perform alignment before, after or both, before and after aggregation.
#' @param alignMethod           Character, Alignment method (see \code{MALDIquant::alignSpectra}.)
#'
#' @return List of preprocessed \code{MALDIquant::MassSpectrum}
#' @export
#' @examples 
#' test_spectra_proc <- AssaySupport::preprocess_spectra(spec = test_spectra,
#'                                                       smooth_meth = "SavitzkyGolay",
#'                                                       smooth_halfWindowSize  = 5,
#'                                                       norm_meth = "TIC",
#'                                                       filter_spectra = NA,
#'                                                       baseline_meth = "TopHat",
#'                                                       align = "none")
#' 
preprocess_spectra <- function(spec = spectra, 
                               smooth_meth = "SavitzkyGolay",
                               smooth_halfWindowSize  = 20,
                               norm_meth = "TIC",
                               filter_spectra = "Cal|D1|D2|CAL",
                               baseline_meth = "TopHat",
                               avg_method = "mean",
                               align = c("none", "before", "after", "both"),
                               alignMethod = "linear") {
  align <- match.arg(align)
  spec_names <- names(spec)
  cat("\n", AlzTools::timeNowHM(), "Preprocessing spectra...\n")
  
  spec <- MALDIquant::calibrateIntensity(spec, method = norm_meth)
  spec <- MALDIquant::smoothIntensity(spec, method = smooth_meth, halfWindowSize  = smooth_halfWindowSize)
  names(spec) <- spec_names
  
  if(any(c("before", "both") %in% align)) {
    aligned_specs <- list()
   for(name in unique(spec_names)) { 
     current_specset <- spec[which(names(spec) == name)]
     current_ref <- MALDIquant::referencePeaks(MALDIquant::detectPeaks(current_specset, SNR = 2, method = "MAD"))
     current_specset_aligned <- MALDIquant::alignSpectra(current_specset,
                                                 warpingMethod = alignMethod, 
                                                 reference = current_ref
                                                 )
     aligned_specs <- append(aligned_specs, current_specset_aligned)
   } 
  }
  
  if(!is.na(filter_spectra)){
    spec <- spec[!grepl(filter_spectra, names(spec))]
  }
  if(!is.na(avg_method)) {
    spec <- MALDIquant::averageMassSpectra(spec, labels = spec_names, method = avg_method)
    spec_names <- unique(spec_names)
    names(spec) <- spec_names
  }
  spec <- MALDIquant::removeBaseline(spec, method = baseline_meth)
  names(spec) <- spec_names
  
  if(any(c("after", "both") %in% align)) {
    spec <- MALDIquant::alignSpectra(spec, reference = MALDIquant::referencePeaks(MALDIquant::detectPeaks(spec)))
    names(spec) <- spec_names
  }
  
  return(spec)
}
