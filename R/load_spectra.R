
#' load spectra of the bruker flex format
#'
#' @param Dir         Character, directory containing spectra 
#' @param filter      Character vector containing spectra names to be filtered
#' @param nameSpectra Logical, name spectra according to there parent folder (default)
#'
#' @return list of MALDIquant spectra
#' @export
#' 
#' @importFrom svMisc progress
load_spectra <- function(Dir, filter = c("Cal", "cal", "CAL"), nameSpectra = TRUE) {
  
  
  # get names of all anaylses
  analyses <- basename(list.dirs(Dir, recursive = F))
  
  # filter Calibration spectra
  analyses <- analyses[which(!analyses %in% filter)]
  
  
  
  # determine number of measured spots
  total_n <- vector("numeric", 1)
  for(i in 1:length(analyses)) {
    n <- length(list.dirs(file.path(Dir, analyses[i]), recursive = F))
    total_n <- total_n + n
    
  }
  
  spectra <- vector("list", length = total_n)
  counter <- 0
  cat(timeNowHM(), "Loading spectra...\n\n")
  for(i in analyses) {
    
    path <- file.path(Dir, i)
    spots <- list.files(path, recursive = T)[grepl(pattern = "fid", x = list.files(path, recursive = T))]
    for(spot in spots) {
      counter <- counter + 1
      spec <- MALDIquantForeign::importBrukerFlex(path = file.path(path, spot), verbose = F )
      #metaData(spec[[1]]) <- list(name = i)
      spectra[[counter]] <- spec[[1]]
      
      if(nameSpectra) {
        names(spectra)[counter] <- i
      }
      
      svMisc::progress(counter/total_n*100)
    }
    
  }
  cat("\n")

return(spectra)
}
