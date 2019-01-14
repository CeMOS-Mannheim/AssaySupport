#' Generate result data.frame containing species found in supplied list
#'
#' @param peak_df           Data.frame containing peak data as generate by \code{generate_peakdf}
#' @param speciesdf         Data.frame, containing columns:
#'                            "Substrate"  Character, the substrate from which the species was generated (used as grouping variable )
#'                            "Species"    Character, the name of the species 
#'                            "mz"         Numeric, m/z of species.
#' @param tol               Numeric, tolerance for assignment.
#' @param totalInt_varNames Character vector, names of variables to sum up for total intensity (e.g. part of product name: "Ab")
#' @param separateIdInto    Character vector, names of new columns into which the ID varibales (e.g. sample name from instrument seperated by \code{sep})
#' @param sep               Character, seperation character for sample names in "ID" column
#'
#' @return
#' data.frame only containing peak information which could be annotated by species given in \code{speciesdf}
#' @export
#'
#' @examples 
#' head(generate_resultdf(test_peak_df, test_Ablist, tol = 5))

generate_resultdf <- function(peak_df, 
                              speciesdf, 
                              tol, 
                              totalInt_varNames = c("Erb", "Neu", "Ab"), 
                              separateIdInto = c("Substrate", "Type", "Condition"),
                              sep = "_") {
  cat("\n", AlzTools::timeNowHM(), "Generating result data.frame...\n")
  label_df <- peak_df %>%
    tidyr::separate(ID, 
                    into = separateIdInto, 
                    sep = sep) %>%
    AssaySupport::assign_species(peak_df = ., 
                                 speciesdf = speciesdf, 
                                 tol = tol, 
                                 mzcol = "mz")
  
  res_df <- label_df %>%
    dplyr::filter(!is.na(species)) %>%
    dplyr::group_by(Substrate, Condition, Type) %>%
    dplyr::mutate(total.Ab = sum(int[grepl(paste(totalInt_varNames, collapse = "|"), species)]),
                  norm.tot = int/total.Ab*100)
  
  return(res_df)
}
