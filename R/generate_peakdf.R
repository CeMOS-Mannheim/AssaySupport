#' generate peak statistic data.frame
#'
#' @param spectra        \code{MALDIquant::MassSpectrum}
#' @param pick_meth      Character, peak picking method
#' @param SNR            numeric, SNR value
#' @param halfWindowSize Integer, halfWindowSize
#' @param binpeaks       Logical, perform binning
#'
#' @return data.frame containing peak statistics
#' @export
#' 
#' @examples 
#' head(generate_peakdf(test_spectra_proc))

generate_peakdf <- function(spectra, pick_meth = "SuperSmoother", halfWindowSize = 10, SNR = 3, binpeaks = FALSE) {
  
  
  
  peaks <- MALDIquant::detectPeaks(spectra, method = pick_meth, SNR = SNR, halfWindowSize = halfWindowSize)
  if(MALDIquant::isMassSpectrum(spectra)) {
    spectra <- list(spectra)
    peaks <- list(peaks)
  }
  
  if(is.null(names(spectra))) {
    names(spectra) <- "noName"
  }
  
  if(binpeaks) {
    peaks <- MALDIquant::binPeaks(peaks)
  }
  
  names(peaks) <- names(spectra)
  
  peak_df <- data.frame(ID = character(),
                        plotIdx <- integer(),
                        mz = double(),
                        int = double(),
                        SNR = double())
  
  for(i in 1:length(peaks)) {
    df <- data.frame(ID = paste0(names(peaks[i])),
                     plotIdx = i,
                     mz = MALDIquant::mass(peaks[[i]]),
                     int = MALDIquant::intensity(peaks[[i]]),
                     SNR = MALDIquant::snr(peaks[[i]]))
    peak_df <- rbind(df, peak_df)
  }
  return(peak_df)
}