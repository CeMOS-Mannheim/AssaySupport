#' Convert vector of m/z values to a peakdf as used by many function
#'
#' @param mz Numerical, mz values
#' @param ID Character, ID as created by FlexControl in format "SUBSTRATE_TYPE_CONDITION"
#'
#' @return
#' peakdf
#' @export

vector2peakdf <- function(mz, ID = "C99_0_0") {
  df <- data.frame(ID = ID,
                   plotIdx = 1,
                   mz = mz,
                   int = 1,
                   SNR = 3)
  return(df)
}
