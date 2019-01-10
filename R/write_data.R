#' write peak statistics as xlsx file to disk
#'
#' @param labels_df n
#' @param parentDir n
#' @param namePrefix n
#'
#' @export
write_peak_list <-function(labels_df = label_df, 
                           parentDir = parentDir, 
                           namePrefix = namePrefix) {
  labels_df %>% 
    filter(!is.na(species)) %>% 
    as.data.frame() %>%
    xlsx::write.xlsx2(file = file.path(parentDir, paste0(namePrefix, "_", basename(parentDir),"_peak_list.xlsx")))
}

#' write spectra as MSD to disk
#'
#' @param spectra List or single \code{MALDIquant::MassSpectrum}, spectra to export.
#' @param parentDir Character, folder to write spectra
#' @param namePrefix Character, prefix to add to filenames
#' @param spectraNames Character vector, names of spectra by default the \code{spectra} are expected to be named.
#' @param peaks List or single \code{MALDIquant::MassPeaks}, peaks to be added to msd file. Set to NA if no peaks should be added.
#' @param annotation_df data.frame containing assigned species. Must contain columns: "species", "mz.theo", "plotIdx", "int", "mz"
#' 
#' @importFrom magrittr %>%  
#' @importFrom dplyr filter select pull
#'
#' @export
write_msd <- function(spectra = spectra_align, 
                      parentDir = parentDir, 
                      spectraNames = unique(paste0(names(spectra))),
                      namePrefix = namePrefix,
                      peaks = peaks,
                      annotation_df = NA) {
  if(is.data.frame(annotation_df) && is.na(peaks[1])) {
    stop("Annotation only possible if peaks are supplied.")
  }
  cat("\n", AlzTools::timeNowHM(), "Writing avg spectra to disk...\n")
  # export spectra
  for(i in 1:length(spectra_align)) {
    current_file <- paste0(parentDir, namePrefix, "_", spectraNames[i], "_avg.msd")
    if(!is.na(peaks[1])) {
      MALDIquantForeign::exportMsd(x = spectra_align[[i]], 
                                   peaks = peaks[[i]], 
                                   file = current_file, 
                                   force = TRUE)
    } else {
      MALDIquantForeign::exportMsd(x = spectra_align[[i]], 
                                   file = current_file, 
                                   force = TRUE)
    }
    if(is.data.frame(annotation_df) ){
      current_df <- annotation_df %>%
                      dplyr::ungroup() %>%
                      dplyr::filter(plotIdx == i, !is.na(species)) %>%
                      dplyr::select(species, mz.theo, int, mz)
                    
      if(nrow(current_df < 1)) {
        warning("No labels found for file", current_file, "\n")
        next()
      }
      for(j in 1:nrow(current_df))
      annotate_msd_peak(file = current_file, 
                        label = current_df[j, "species"] %>% dplyr::pull(), 
                        peaks = c(current_df[j, "mz",] %>% dplyr::pull(),  current_df[j, "int"] %>% dplyr::pull()), 
                        peakIdx = NA, 
                        theoMz =  current_df[j, "mz.theo"] %>% dplyr::pull())
      
    }
    
  }
  
}
