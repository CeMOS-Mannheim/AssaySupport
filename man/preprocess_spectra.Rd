% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{preprocess_spectra}
\alias{preprocess_spectra}
\title{Preprocess spectra}
\usage{
preprocess_spectra(
  data,
  smooth_meth = "SavitzkyGolay",
  smooth_halfWindowSize = 20,
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
  allowNoMatch = TRUE
)
}
\arguments{
\item{data}{List of \code{MALDIquant::MassSpectrum}}

\item{smooth_meth}{Character, smoothing method see \code{MALDIquant::smoothIntensity()}}

\item{smooth_halfWindowSize}{Numeric, halfWindowSize for smoothing method see \code{MALDIquant::smoothIntensity()}}

\item{norm_meth}{Character, normalization method see \code{MALDIquant::calibrateIntensity()}.
Additional to the options that MALDIquant provides for normalization it is also possible to normalize to a internal standard ("IS").
To use this option \code{IS_mass} and \code{IS_tol} has to be set .}

\item{lockMass}{Numeric, mass for lock mass recalibration. See \code{generate_resultdf}. 
Set to \code{NA} if no mass recalibartion should be performed.}

\item{lockMass_tol}{Numeric, tolerance for lock mass species assignment.}

\item{IS_mass}{Character, name of the species to use as internal standard for normalization. See \code{generate_resultdf}.}

\item{IS_tol}{Numeric, tolerance for internal standard species assignment.}

\item{filter_spectra}{Character vector, regex of patterns to exclude spectra (like calibration spectra in same folder)}

\item{baseline_meth}{Character, baseline removal method see \code{MALDIquant::removeBaseline}}

\item{avg_method}{Character, aggregation method used to generate average spectra. See \code{MALDIquant::averageMassSpectra}.}

\item{align}{Chacter vector, perform alignment before, after or both, before and after aggregation.}

\item{align_Method}{Character, Alignment method (see \code{MALDIquant::alignSpectra}).}

\item{align_SNR}{Numeric, SNR for peak picking (see \code{MALDIquant::alignSpectra}).}

\item{ISlockMass_SNR}{Numeric, SNR for detection of the IS/lock mass.}

\item{tolppm}{Logical, all tolerances in ppm instead of Dalton.}

\item{align_pickMeth}{Character, Peak picking method (see \code{MALDIquant::alignSpectra}).}

\item{align_minFreq}{Character, minimal peak frequency (see \code{MALDIquant::alignSpectra}).}

\item{align_tol}{Character, tolerance to consider peak as identical (see \code{MALDIquant::alignSpectra}).}

\item{allowNoMatch}{Logical, spectra with no match for IS normalization or lock mass calibration will be excluded. If \code{FALSE} method will stop if no match can be found.}
}
\value{
List of preprocessed \code{MALDIquant::MassSpectrum}
}
\description{
Preprocess spectra
}
\examples{
test_spectra_proc <- AssaySupport::preprocess_spectra(data = test_spectra,
                                                      smooth_meth = "SavitzkyGolay",
                                                      smooth_halfWindowSize  = 5,
                                                      norm_meth = "TIC",
                                                      filter_spectra = NA,
                                                      baseline_meth = "TopHat",
                                                      align = "none")

}
