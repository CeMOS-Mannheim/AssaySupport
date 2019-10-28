#' Preprocess spectra
#'
#' @param spec                  List of \code{MALDIquant::MassSpectrum}
#' @param smooth_meth           Character, smoothing method see \code{MALDIquant::smoothIntensity()}
#' @param smooth_halfWindowSize Numeric, halfWindowSize for smoothing method see \code{MALDIquant::smoothIntensity()}
#' @param norm_meth             Character, normalization method see \code{MALDIquant::calibrateIntensity()}.
#'                              Additional to the options that MALDIquant provides for normalization it is also possible to normalize to a internal standard ("IS").
#'                              To use this option \code{IS_species} and \code{IS_tol} has to be set and a \code{species_df} has to be provided (see \code{AssaySupport::generate_assigndf()}).
#' @param lockMass_species      Character, name of the species to use as lock mass. See \code{generate_resultdf}. 
#'                              Set to \code{NA} if no mass recalibartion should be performed.
#' @param lockMass_tol          Numeric, tolerance for lock mass species assignment.
#' @param lockMass_tolPPM       Logical, use the \code{lockMass_tol} as ppm.
#' @param IS_species            Character, name of the species to use as internal standard for normalization. See \code{generate_resultdf}.
#' @param IS_tol                Numeric, tolerance for internal standard species assignment.
#' @param IS_tolPPM             Logical, use the \code{IS_tol} as ppm.
#' @param filter_spectra        Character vector, regex of patterns to exclude spectra (like calibration spectra in same folder)
#' @param baseline_meth         Character, baseline removal method see \code{MALDIquant::removeBaseline}
#' @param avg_method            Character, aggregation method used to generate average spectra. See \code{MALDIquant::averageMassSpectra}.
#' @param align                 Chacter vector, perform alignment before, after or both, before and after aggregation.
#' @param align_Method          Character, Alignment method (see \code{MALDIquant::alignSpectra}).
#' @param align_SNR             Numeric, SNR for peak picking (see \code{MALDIquant::alignSpectra}).
#' @param align_pickMeth        Character, Peak picking method (see \code{MALDIquant::alignSpectra}).
#' @param align_minFreq         Character, minimal peak frequency (see \code{MALDIquant::alignSpectra}).
#' @param align_tol             Character, tolerance to consider peak as identical (see \code{MALDIquant::alignSpectra}).
#' @param species_df            Data.frame, see \code{generate_assigndf}.
#' @param allowNoMatch          Logical, spectra with no match for IS normalization or lock mass calibration will be excluded. If \code{FALSE} method will stop if no match can be found.
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
                               norm_meth = c("TIC", "IS", "none"),
                               lockMass_species = NA,
                               lockMass_tol = 250,
                               lockMass_tolPPM = TRUE,
                               IS_species = lockMass_species,
                               IS_tol = lockMass_tol,
                               IS_tolPPM = lockMass_tolPPM,
                               filter_spectra = "Cal|D1|D2|CAL",
                               baseline_meth = "TopHat",
                               avg_method = c("mean", "median", "sum", NA),
                               align = c("none", "before", "after", "both"),
                               align_Method = "linear",
                               align_SNR = 2,
                               align_pickMeth = "SuperSmoother",
                               align_minFreq = 0.25,
                               align_tol = 0.01,
                               species_df = C99T_commonSpecies,
                               allowNoMatch = TRUE) {
  norm_meth <- match.arg(norm_meth)
  avg_method = match.arg(avg_method)
  align <- match.arg(align)
  
  if(norm_meth == "IS") {
    if(is.na(IS_species) | dim(species_df)[2] < 1) {
      stop("IS_species and species_df must be provided for IS normalization.\n")
    }
    if(!(IS_species %in% species_df$Species)) {
      stop("IS_species not containied in species_df.\n")
    }
  }
  
  
  spec_names <- names(spec)
  
  cat("\n", AlzTools::timeNowHM(), "Preprocessing spectra...\n")
  
  spec <- MALDIquant::smoothIntensity(spec, method = smooth_meth, halfWindowSize  = smooth_halfWindowSize)
  spec <- MALDIquant::removeBaseline(spec, method = baseline_meth)
  names(spec) <- spec_names
  
  switch(norm_meth,
         TIC = {
           spec <- MALDIquant::calibrateIntensity(spec, 
                                                  method = norm_meth)
         }, 
         median = {
           spec <- MALDIquant::calibrateIntensity(spec, 
                                                  method = norm_meth)
         },
         PQM = {
           spec <- MALDIquant::calibrateIntensity(spec, 
                                                  method = norm_meth)
         },
         IS = {
           peak_df <- generate_peakdf(spectra = spec, 
                                      SNR = 3)
           norm_fac <- getNormFactors(peaksdf = peak_df,
                                      speciesdf = species_df, 
                                      targetSpecies = IS_species, 
                                      tol = IS_tol, 
                                      tolppm = IS_tolPPM, 
                                      allowNoMatch = allowNoMatch)
           
           spec <- normalizeByFactor(spec[norm_fac$specIdx], norm_fac$norm_factor)
           spec_names <- spec_names[norm_fac$specIdx]
         }, 
         none = {
           # don't do anything
         })
  names(spec) <- spec_names
  
  if(any(c("before", "both") %in% align)) {
    aligned_specs <- list()
    for(name in unique(spec_names)) { 
      current_specset <- spec[which(names(spec) == name)]
      
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
    spec <- aligned_specs
  }
  
  if(!is.na(lockMass_species)) {
    peak_df <- generate_peakdf(spectra = spec, 
                               SNR = 3)
    mzshift <- getMzShift(peak_df, 
                          speciesdf = species_df, 
                          targetSpecies = lockMass_species, 
                          tol = lockMass_tol, 
                          tolppm = lockMass_tolPPM,
                          allowNoMatch = allowNoMatch)
    spec <- shiftMassAxis(spec[mzshift$specIdx], mzshift$mzshift)
    spec_names <- spec_names[mzshift$specIdx]
  }
  
  names(spec) <- spec_names
  
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
    spec <- MALDIquant::alignSpectra(spec,
                                     warpingMethod = align_Method, 
                                     reference = MALDIquant::referencePeaks(
                                       MALDIquant::detectPeaks(spec, 
                                                               SNR = align_SNR, 
                                                               method = align_pickMeth),
                                       minFrequency = align_minFreq, 
                                       tolerance = align_tol))
    names(spec) <- spec_names
  }
  
  return(spec)
}
