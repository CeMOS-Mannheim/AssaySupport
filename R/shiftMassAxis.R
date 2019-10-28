#' Shift mass axis/Linear point correction
#'
#' @param spec    MALDIquant::MassSpectrum or MALDIquant::MassPeaks
#' @param mzdiff  numeric, correction value for mass axis.
#'
#' @return
#' MALDIquant::MassSpectrum or MALDIquant::MassPeaks with shifted mass axis
#' @export
shiftMassAxis <- function(spec, mzdiff) {
  if(isMassSpectrum(spec) || isMassPeaks(spec)) {
    spec@mass <- spec@mass + mzdiff
    return(spec)  
  }
  
  if(isMassSpectrumList(spec) || isMassPeaksList(spec)) {
    if(!(length(spec) == length(mzdiff))) {
      stop("length(spec) != length(mzdiff) !\n")
    }
    for(i in 1:length(spec)) {
      spec[[i]]@mass <- spec[[i]]@mass + mzdiff[i]
    }
    return(spec)
  }
  stop("spec needs to be a MALDIquant::MassSpectrum or MALDIquant::MassPeaks or a list of these. \n")
}

