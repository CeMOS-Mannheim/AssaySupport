#' generate peak statistic data.frame
#'
#' @param data           \code{MALDIquant::MassSpectrum} or \code{MALDIquant::MassPeaks}
#' @param pick_meth      Character, peak picking method. Only applies is no peaks are supplied.
#' @param SNR            numeric, SNR value. Only applies is no peaks are supplied.
#' @param halfWindowSize Integer, halfWindowSize. Only applies is no peaks are supplied.
#' @param binpeaks       Logical, perform binning.
#'
#' @return data.frame containing peak statistics
#' @export
#' 
#' @examples 
#' head(generate_peakdf(test_spectra_proc))

generate_peakdf <- function(data, 
                            pick_meth = "SuperSmoother",
                            halfWindowSize = 10, 
                            SNR = 3, 
                            binpeaks = FALSE,
                            binTol = 50e-6) {
  if(MALDIquant::isMassPeaksList(data) || MALDIquant::isMassPeaks(data)) {
    peaks <- data
  } else {
    peaks <- MALDIquant::detectPeaks(data, method = pick_meth, SNR = SNR, halfWindowSize = halfWindowSize)
  }
  
  if(MALDIquant::isMassSpectrum(data) || MALDIquant::isMassPeaks(data)) {
    data <- list(data)
    peaks <- list(peaks)
  }
  
  if(is.null(names(data))) {
    names(data) <- "noName"
  }
  
  if(binpeaks) {
    peaks <- MALDIquant::binPeaks(peaks, tolerance = binTol)
  }
  
  names(peaks) <- names(data)
  
  peak_df <- data.frame(ID = character(),
                        plotIdx = integer(),
                        peakIdx = character(),
                        mz = double(),
                        int = double(),
                        SNR = double())
  
  for(i in 1:length(peaks)) {
    df <- data.frame(ID = paste0(names(peaks[i])),
                     plotIdx = i,
                     peakIdx = paste0(i, "_", 1:length(MALDIquant::mass(peaks[[i]]))),
                     mz = MALDIquant::mass(peaks[[i]]),
                     int = MALDIquant::intensity(peaks[[i]]),
                     SNR = MALDIquant::snr(peaks[[i]]))
    peak_df <- rbind(df, peak_df)
  }
  return(peak_df)
}
