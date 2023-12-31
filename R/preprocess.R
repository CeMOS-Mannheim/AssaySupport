#' Preprocess spectra
#'
#' @param data                  List of \code{MALDIquant::MassSpectrum}
#' @param smooth_meth           Character, smoothing method see \code{MALDIquant::smoothIntensity()}
#' @param smooth_halfWindowSize Numeric, halfWindowSize for smoothing method see \code{MALDIquant::smoothIntensity()}
#' @param norm_meth             Character, normalization method see \code{MALDIquant::calibrateIntensity()}.
#'                              Additional to the options that MALDIquant provides for normalization it is also possible to normalize to a internal standard ("IS").
#'                              To use this option \code{IS_mass} and \code{IS_tol} has to be set .
#' @param lockMass              Numeric, mass for lock mass recalibration. See \code{generate_resultdf}. 
#'                              Set to \code{NA} if no mass recalibartion should be performed.
#' @param lockMass_tol          Numeric, tolerance for lock mass species assignment.
#' @param IS_mass            Character, name of the species to use as internal standard for normalization. See \code{generate_resultdf}.
#' @param IS_tol                Numeric, tolerance for internal standard species assignment.
#' @param filter_spectra        Character vector, regex of patterns to exclude spectra (like calibration spectra in same folder)
#' @param baseline_meth         Character, baseline removal method see \code{MALDIquant::removeBaseline}
#' @param avg_method            Character, aggregation method used to generate average spectra. See \code{MALDIquant::averageMassSpectra}.
#' @param align                 Chacter vector, perform alignment before, after or both, before and after aggregation.
#' @param align_Method          Character, Alignment method (see \code{MALDIquant::alignSpectra}).
#' @param align_SNR             Numeric, SNR for peak picking (see \code{MALDIquant::alignSpectra}).
#' @param ISlockMass_SNR        Numeric, SNR for detection of the IS/lock mass.
#' @param tolppm                Logical, all tolerances in ppm instead of Dalton.
#' @param align_pickMeth        Character, Peak picking method (see \code{MALDIquant::alignSpectra}).
#' @param align_minFreq         Character, minimal peak frequency (see \code{MALDIquant::alignSpectra}).
#' @param align_tol             Character, tolerance to consider peak as identical (see \code{MALDIquant::alignSpectra}).
#' @param allowNoMatch          Logical, spectra with no match for IS normalization or lock mass calibration will be excluded. If \code{FALSE} method will stop if no match can be found.
#'
#' @return List of preprocessed \code{MALDIquant::MassSpectrum}
#' @export
#' @examples 
#' test_spectra_proc <- AssaySupport::preprocess_spectra(data = test_spectra,
#'                                                       smooth_meth = "SavitzkyGolay",
#'                                                       smooth_halfWindowSize  = 5,
#'                                                       norm_meth = "TIC",
#'                                                       filter_spectra = NA,
#'                                                       baseline_meth = "TopHat",
#'                                                       align = "none")
#' 
preprocess_spectra <- function(data, 
                               smooth_meth = "SavitzkyGolay",
                               smooth_halfWindowSize  = 20,
                               norm_meth = c("TIC", "IS", "none"),
                               lockMass = NA,
                               lockMass_tol = 250,
                               IS_mass = lockMass,
                               IS_tol = lockMass_tol,
                               filter_spectra = NA,
                               baseline_meth = "TopHat",
                               avg_method = c("mean", "median", "sum", NA),
                               align = c("none", "before", "after", "both"),
                               align_Method = "linear",
                               align_SNR = 2,
                               ISlockMass_SNR = 3, 
                               tolppm = TRUE,
                               align_pickMeth = "SuperSmoother",
                               align_minFreq = 0.25,
                               align_tol = 0.01,
                               allowNoMatch = TRUE) {
  norm_meth <- match.arg(norm_meth)
  avg_method = match.arg(avg_method)
  align <- match.arg(align)
  
  if(norm_meth == "IS") {
    if(is.na(IS_mass)) {
      stop("IS_mass must be provided for IS normalization.\n")
    }
  }
  
  
  spec_names <- names(data)
  
  cat("\n", timeNowHM(), "Preprocessing spectra...\n")
  
  data <- MALDIquant::smoothIntensity(data, method = smooth_meth, halfWindowSize  = smooth_halfWindowSize)
  data <- MALDIquant::removeBaseline(data, method = baseline_meth)
  names(data) <- spec_names
  
  switch(norm_meth,
         TIC = {
           data <- MALDIquant::calibrateIntensity(data, 
                                                  method = norm_meth)
         }, 
         median = {
           data <- MALDIquant::calibrateIntensity(data, 
                                                  method = norm_meth)
         },
         PQM = {
           data <- MALDIquant::calibrateIntensity(data, 
                                                  method = norm_meth)
         },
         IS = {
           peak_df <- generate_peakdf(data, 
                                      SNR = ISlockMass_SNR)
           norm_fac <- getNormFactors(peaksdf = peak_df,
                                      speciesdf = NA, 
                                      targetSpecies = IS_mass, 
                                      tol = IS_tol, 
                                      tolppm = tolppm, 
                                      allowNoMatch = allowNoMatch)
           
           data <- normalizeByFactor(data[norm_fac$specIdx], norm_fac$norm_factor)
           spec_names <- spec_names[norm_fac$specIdx]
         }, 
         none = {
           # don't do anything
         })
  names(data) <- spec_names
  
  if(any(c("before", "both") %in% align)) {
    aligned_specs <- list()
    for(name in unique(spec_names)) { 
      current_specset <- data[which(names(data) == name)]
      
      current_ref <- MALDIquant::referencePeaks(
                                MALDIquant::detectPeaks(current_specset, 
                                                        SNR = align_SNR, 
                                                        method = align_pickMeth),
                                minFrequency = align_minFreq, tolerance = align_tol)
      current_specset_aligned <- MALDIquant::alignSpectra(current_specset,
                                                          warpingMethod = align_Method, 
                                                          reference = current_ref)
      aligned_specs <- append(aligned_specs, current_specset_aligned)
    } 
    data <- aligned_specs
  }
  
  if(!is.na(lockMass)) {
    peak_df <- generate_peakdf(data, 
                               SNR = ISlockMass_SNR)
    mzshift <- getMzShift(peak_df, 
                          speciesdf = NA, 
                          targetSpecies = lockMass, 
                          tol = lockMass_tol, 
                          tolppm = tolppm,
                          allowNoMatch = allowNoMatch)
    data <- shiftMassAxis(data[mzshift$specIdx], mzshift$mzshift)
    spec_names <- spec_names[mzshift$specIdx]
  }
  
  names(data) <- spec_names
  
  if(!is.na(filter_spectra)){
    data <- data[!grepl(filter_spectra, names(data))]
  }
  if(!is.na(avg_method)) {
    data <- MALDIquant::averageMassSpectra(data, labels = spec_names, method = avg_method)
    spec_names <- unique(spec_names)
    names(data) <- spec_names
  }
  data <- MALDIquant::removeBaseline(data, method = baseline_meth)
  names(data) <- spec_names
  
  if(any(c("after", "both") %in% align)) {
    data <- MALDIquant::alignSpectra(data,
                                     warpingMethod = align_Method, 
                                     reference = MALDIquant::referencePeaks(
                                       MALDIquant::detectPeaks(data, 
                                                               SNR = align_SNR, 
                                                               method = align_pickMeth),
                                       minFrequency = align_minFreq, 
                                       tolerance = align_tol))
    names(data) <- spec_names
  }
  
  return(data)
}
