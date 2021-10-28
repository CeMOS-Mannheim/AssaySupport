#' Generate result data.frame containing species found in supplied list
#'
#' @param peak_df           Data.frame containing peak data as generate by \code{generate_peakdf}
#' @param speciesdf         Data.frame, containing columns:
#'                            "Substrate"  Character, the substrate from which the species was generated (used as grouping variable )
#'                            "Species"    Character, the name of the species 
#'                            "mz"         Numeric, m/z of species.
#' @param tol               Numeric, tolerance for assignment.
#' @param highMz      Numeric, if provided a second tolerance (\code{highMzTol}) for \code{mz >= highMz} will be used.
#' @param highMzTol   Numeric, tolerance for assignment for \code{mz >= highMz} either in Da or ppm.
#' @param tolppm            Logical, use relative tolerance in ppm. Supply \code{tol} as e.g. 100e-6
#' @param nHits             Numeric, maximal number of hits whithin tolerance. Assignment only, no quantitative analysis.
#' @param totalInt_varNames Character vector, names of variables to sum up for total intensity (e.g. part of product name: "Ab")
#' @param separateIdInto    Character vector, names of new columns into which the ID varibales (e.g. sample name from instrument seperated by \code{sep})
#' @param subsetcol         Character, column of \code{peak_df} used for subsetting (usually the substrate used in reaction.)
#' @param sep               Character, seperation character for sample names in "ID" column
#'
#' @return
#' data.frame only containing peak information which could be annotated by species given in \code{speciesdf}
#' @export
#' 
#' @importFrom tibble tibble
#'
#' @examples 
#' head(generate_resultdf(test_peak_df, test_Ablist, tol = 5))

generate_resultdf <- function(peak_df, 
                              speciesdf, 
                              tol,
                              highMz = NA,
                              highMzTol = NA,
                              tolppm = TRUE,
                              nHits = 1,
                              totalInt_varNames = c("Erb", "Neu", "Ab"), 
                              separateIdInto = c("Substrate", "Type", "Condition"),
                              subsetcol =  separateIdInto[1],
                              sep = "_") {
  if(missing(tol)) {
    stop("No tolerance provided. Please set 'tol'.\n")
  }
  if(missing(peak_df)) {
    stop("No peak_df provided.\n")
  }
  if(missing(speciesdf)) {
    stop("No speciesdf provided.\n")
  }
  if(!subsetcol %in% colnames(speciesdf)) {
    stop("subsetcol not found in colnames of speciesdf.\n")
  }
  
  cat("\n", timeNowHM(), "Generating result data.frame...\n")
  
  
  
  res_df <- tibble()
  peak_df <- peak_df %>%
    tidyr::separate(ID, 
                    into = separateIdInto, 
                    sep = sep)
  
  if(!all(peak_df[,subsetcol] %in% speciesdf[,subsetcol])) {
    stop("Not all entries in peak_df[,subsetcol] were found in speciesdf[,subsetcol].\n",
         "Entries in speciesdf: ", paste(unique(speciesdf[,subsetcol]), collapse = ", "), ".\n",
         "Entries in peak_df: ", paste(unique(peak_df[,subsetcol]), collapse = ", "), ".")
  }
  
  Idx <- peak_df %>%
    pull(plotIdx) %>%
    unique()
  
  if(nHits == 1) {
    label_df <- peak_df %>%
      AssaySupport::assign_species(peak_df = ., 
                                   subsetcol = subsetcol,
                                   speciesdf = speciesdf, 
                                   tol = tol, 
                                   highMz = highMz,
                                   highMzTol = highMzTol,
                                   tolppm = tolppm,
                                   mzcol = "mz")
  } else {
    speciesdf_orig <- speciesdf
    for(i in 1:length(Idx)) {
      speciesdf <- speciesdf_orig
      current_res_df <- tibble()
      for(j in 1:nHits) {
        label_df <- peak_df %>%
          filter(plotIdx == Idx[i]) %>%
          AssaySupport::assign_species(peak_df = ., 
                                       subsetcol = subsetcol,
                                       speciesdf = speciesdf, 
                                       tol = tol, 
                                       highMz = highMz,
                                       highMzTol = highMzTol,
                                       tolppm = tolppm,
                                       mzcol = "mz")
        
        current_res_df <- bind_rows(current_res_df, label_df)
      speciesdf <- speciesdf %>%
        filter(!Species %in% current_res_df$species)
      }
      res_df <- bind_rows(res_df, current_res_df)
    }
    
  }
  
  if(nHits > 1) {
    return(res_df)
  }
  
  res_df <- label_df %>%
    dplyr::filter(!is.na(species)) %>%
    dplyr::group_by(plotIdx) %>%
    dplyr::mutate(total.Ab = sum(int[grepl(paste(totalInt_varNames, collapse = "|"), species)]),
                  norm.tot = int/total.Ab*100)
  
  return(res_df)
}
