#' Read csv of overview spectra from SCiLS
#'
#' @param path Character, path to folder with .csv files
#' @param type Character, type of spectrum to import.
#'
#' @return
#' List of MALDIquant::MassSpectrum
#' @importFrom tibble as_tibble
#' @importFrom dplyr filter mutate pull %>%
#' @importFrom tidyr separate
#' @importFrom stringr str_detect str_remove
#' @export

readScilsCsv <- function(path, type = c("mean", "skyline", "variance")) {
  type <- match.arg(type)
  
  files <- list.files(path)[grepl(".csv", list.files(path))]
  
  result <- vector("list", length = length(files))
  for(i in 1:length(files)) {
    spec_csv <- read.csv(file.path(path, files[i]), comment.char = "#")
    
    # prepare the metadata part of the csv
    meta_data <- read.csv(file.path(path, files[i]), nrows = 20, header = FALSE) %>% 
      as_tibble() %>%
      filter(str_detect(V1, "#")) %>%
      separate(V1, into = c("name", "value"), sep = ": ", fill = "right") %>%
      mutate(name = str_remove(name, "# "))
    
    scilsVersion <- meta_data %>% 
      filter(str_detect(name, "Version")) %>%
      separate(name, into = c("rest", "version"), sep = "Version ") %>%
      pull(version)
    
    exportTime <- meta_data %>%
      filter(str_detect(name, "Export time")) %>%
      pull(value)
    
    creationTime <- meta_data %>%
      filter(str_detect(name, "creation time")) %>%
      pull(value)
    
    sourceFile <- meta_data %>%
      filter(str_detect(name, "file")) %>%
      pull(value)
    
    objID <- meta_data %>%
      filter(str_detect(name, "ID")) %>%
      pull(value)
    
    norm <- meta_data %>%
      filter(str_detect(name, "Normalization")) %>%
      pull(value)
    
    type_list <- vector("list", length(type))
    for(j in 1:length(type)) {
      typename <- switch (type[j],
        mean = "intensities",
        skyline = "skyline",
        variance = "variances"
      )
      specname <- paste0(type[j], "_", tools::file_path_sans_ext(files[i]))
      
      type_list[[j]] <- MALDIquant::createMassSpectrum(mass = spec_csv$m.z, 
                                                       intensity = spec_csv[, typename],
                                                       metaData = list(file = sourceFile, 
                                                                       name = specname,
                                                                       exportTime = exportTime,
                                                                       creationTime = creationTime,
                                                                       normalization = norm))
      names(type_list)[j] <- specname
    }
    result[[i]] <- type_list
  }
  return(unlist(result, recursive = FALSE))
}






