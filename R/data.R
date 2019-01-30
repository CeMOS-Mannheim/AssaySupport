#' Sequence of C99 with added 3xFLAG
#' 
#' 131 amino acids. Old sequence as used for reaction in Szaruga et. al., 2017, Cell. 
#' Contains two amino acid long linker region (GT)
"C99_3xFLAG"

#' Sequence of C99 with added 3xFLAG and shorter linker region
#' 
#' 130 amino acids.
#' Removed T from linker region (see \code{C99_3xFLAG}) to shift masses away from contaminations at mass of AICD50-99 and Abeta1-38. 
"C99_deltaT"

#' Sequence of Neurexin
#' 
#' 121 amino acids.
#' Alternative Substrate of GSEC.
"Neurexin"

#' Sequence of ErbB4
#' 
#' 107 amino acids.
#' Alternative Substrate of GSEC.
"ErbB4"

#' Sequence of EphA4
#' 
#' 86 amino acids.
#' Alternative Substrate of GSEC.
"EphA4"

#' Alzheimer's causing mutations in APP (C99).
#' 
#' Collection of APP mutations as used by Matthias Koch (Lucia Chavez Guiterrez Lab, KU Leuven).
#' 
#' @format A data frame with with APP mutations. 4 variables:
#' \describe{
#'   \item{name}{Trivial name of mutation. Mostly location were initial case occured}
#'   \item{code}{Muatation code. First letter is amino acid in WT variant, position of muatation, amino acid in mutated variant}
#'   \item{pos}{Position of muatation. Same as in \code{code}.}
#'   \item{substitute}{amino acid in mutated variant. Same as in \code{code}.}
#'   }
"C99_mutations"

#' Mass list of GSEC reaction related products from different substrates.
#' 
#' @format A data frame with with 4 variables:
#' \describe{
#'   \item{Substrate}{Name of substrate}
#'   \item{Species}{Name of product species, eg. Ab1_38 or AICD49_50}
#'   \item{mz}{Mass of species}
#'   \item{comment}{Free text fild.}
#'   }
"test_Ablist"

#' Spectra (\code{MALDIquant::MassSpectrum}) of GSEC reactions as examples.
"test_spectra"

#' Preprocessed spectra (\code{MALDIquant::MassSpectrum}) of GSEC reactions as examples
"test_spectra_proc"

#' Peaks of (\code{test_spectra_proc}) as \code{MALDIquant::MassPeaks}
"test_peaks"

#' Example peak list of \code{test_spectra}
#' 
#' @format A data frame with with 5 variables:
#' \describe{
#'   \item{ID}{Sample name as defined at istrument when data was recorded.}
#'   \item{plotIdx}{Index of spectrum.}
#'   \item{mz}{Mass of peak.}
#'   \item{int}{Intensity of peak.}
#'   \item{SNR}{Signal-to-Noise Ratio of peak.}
#'   }
"test_peak_df"


