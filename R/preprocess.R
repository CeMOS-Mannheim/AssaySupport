#' Preprocess spectra
#'
#' @param spectra spectra
#' @param smooth_meth smoothing method
#' @param smooth_halfWindowSize n
#' @param norm_meth n
#' @param filter_spectra n 
#' @param baseline_meth n
#'
#' @return preprocessed spectra
#' @export
preprocess_spectra <- function(spec = spectra, 
                               smooth_meth = "SavitzkyGolay",
                               smooth_halfWindowSize  = 20,
                               norm_meth = "TIC",
                               filter_spectra = "Cal|D1|D2|CAL",
                               baseline_meth = "TopHat",
                               align = FALSE) {
  cat("\n", AlzTools::timeNowHM(), "Preprocessing spectra...\n")
  spectra_proc <- MALDIquant::smoothIntensity(spec, method = smooth_meth, halfWindowSize  = smooth_halfWindowSize)
  spectra_proc <- MALDIquant::calibrateIntensity(spectra_proc, method = norm_meth)
  names(spectra_proc) <- paste0(names(spectra))
  spectra_proc_nocal<- spectra_proc[!grepl(filter_spectra,names(spectra_proc))]
  spectra_avg <- MALDIquant::averageMassSpectra(spectra_proc_nocal, labels = names(spectra_proc_nocal), method = "mean")
  spectra_avg <- MALDIquant::removeBaseline(spectra_avg, method = baseline_meth)
  names(spectra_avg) <- unique(names(spectra_proc_nocal))
  if(align) {
    spectra_align <- MALDIquant::alignSpectra(spectra_avg, reference = MALDIquant::referencePeaks(MALDIquant::detectPeaks(spectra_avg)))
    names(spectra_align) <- unique(names(spectra_avg))
    return(spectra_align)
  }
  return(spectra_avg)
}