#' write peak statistics as xlsx file to disk
#'
#' @param labels_df 
#' @param parentDir 
#' @param namePrefix 
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
#' @param spectra 
#' @param parentDir 
#' @param namePrefix 
#'
#' @export
write_msd <- function(spectra = spectra_align, 
                      parentDir = parentDir, 
                      namePrefix = namePrefix,
                      peaks = peaks) {
  cat("\n", AlzTools::timeNowHM(), "Writing avg spectra to disk...\n")
  # export spectra
  for(i in 1:length(spectra_align)) {
    nm <- unique(paste0(names(spectra)))
    if(!is.na(peaks[1])) {
      MALDIquantForeign::exportMsd(x = spectra_align[[i]], 
                                   peaks = peaks[[i]], 
                                   file = paste0(parentDir, namePrefix, "_", nm[i], "_avg.msd"), 
                                   force = TRUE)
    } else {
      MALDIquantForeign::exportMsd(x = spectra_align[[i]], 
                                   file = paste0(parentDir, namePrefix, "_", nm[i], "_avg.msd"), 
                                   force = TRUE)
    }
  }
}