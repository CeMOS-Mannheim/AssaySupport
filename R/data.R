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

#' Common products of GSEC activity on C99_deltaT
#' 
#' Includes all Abetas starting from 1 or 11 and ending in the range of 37:49 as well as AICDs and Species the latter each as [M+H]+ and [M+H2]+2
"C99T_commonSpecies"

#' Products of GSEC activity on C99 including modified formes
#' 
#' Includes all Abetas starting from 1:5, 11, 17 or 19 and ending in the range of 23:49 as well as AICDs and Species the latter each as [M+H]+ and [M+H2]+2 as well as modifications.
"C99_full"

#' Sequence of Neurexin
#' 
#' 121 amino acids.
#' Alternative Substrate of GSEC.
"Neurexin"

#' Sequence of ErbB4
#' 
#' 109 amino acids.
#' Alternative Substrate of GSEC.
"ErbB4"

#' Sequence of EphA4
#' 
#' 86 amino acids.
#' Alternative Substrate of GSEC.
"EphA4"

#' Sequence of N-Cadherin
#' 
#' 109 amino acids.
#' Alternative Substrate of GSEC.
"NCadherin"

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
#' Example of \code{result_df} containing annotations
"test_result_df"

#' Calibration spectrum containing average masses of:
#' name				                avg		    mono
#' ACTH_clip(18-39)		        2466.68		
#' ACTH_clip(1-17)		        2094.42				
#' Somatostatin			          3149.574
#' Ab1-28				              3263.472
#' Ab1-38				              4131.548	4132.5557
#' Ab1-39 				            4230.680	4231.6870
#' Ab1-40 				            4329.811	4330.8182
#' Ab1-42 				            4514.047	4515.0541
#' Ab1-43 				            4615.151	4616.1582
#' Ab1-45 				            4827.440	4828.4473
#' Ab1-46 				            4926.571	4927.5785
#' Ab1-48 				            5140.833	5141.8405
#' Insulin_[M+H]+_avg 		    5734.520
#' Cytochrome_C_[M+2H]2+_avg 	6181.050
#' Myoglobin_[M+2H]2+_avg 		8476.660
#' Ubiquitin_I_[M+H]+_avg 		8565.760
"PAP_cal"

