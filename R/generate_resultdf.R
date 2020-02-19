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
#' @examples 
#' head(generate_resultdf(test_peak_df, test_Ablist, tol = 5))

# generate_resultdf <- function(peak_df, 
#                               speciesdf, 
#                               tol,
#                               highMz = NA,
#                               highMzTol = NA,
#                               tolppm = TRUE,
#                               nHits = 1,
#                               totalInt_varNames = c("Erb", "Neu", "Ab"), 
#                               separateIdInto = c("Substrate", "Type", "Condition"),
#                               subsetcol =  separateIdInto[1],
#                               sep = "_") {
#   cat("\n", AlzTools::timeNowHM(), "Generating result data.frame...\n")
#   res_df <- tibble()
#   
#   for(i in 1:nHits) {
#     label_df <- peak_df %>%
#       tidyr::separate(ID, 
#                       into = separateIdInto, 
#                       sep = sep) %>%
#       AssaySupport::assign_species(peak_df = ., 
#                                    subsetcol = subsetcol,
#                                    speciesdf = speciesdf, 
#                                    tol = tol, 
#                                    highMz = highMz,
#                                    highMzTol = highMzTol,
#                                    tolppm = tolppm,
#                                    mzcol = "mz")
#     if(nHits > 1) {
#       res_df <- bind_rows(res_df, label_df)
#       speciesdf <- speciesdf %>%
#         filter(!Species %in% res_df$species)
#     }
#   }
#   if(nHits > 1) {
#     return(res_df)
#   }
#   
#   res_df <- label_df %>%
#     dplyr::filter(!is.na(species)) %>%
#     dplyr::group_by(plotIdx) %>%
#     dplyr::mutate(total.Ab = sum(int[grepl(paste(totalInt_varNames, collapse = "|"), species)]),
#                   norm.tot = int/total.Ab*100)
#   
#   return(res_df)
# }

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
  cat("\n", AlzTools::timeNowHM(), "Generating result data.frame...\n")
  
  res_df <- tibble()
  peak_df <- peak_df %>%
    tidyr::separate(ID, 
                    into = separateIdInto, 
                    sep = sep)
  
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
        
      res_df <- bind_rows(res_df, label_df)
      speciesdf <- speciesdf %>%
        filter(!Species %in% res_df$species)
      }
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
