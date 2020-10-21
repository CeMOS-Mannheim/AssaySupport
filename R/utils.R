#' get current time in hours and minutes
#'
#' @export

timeNowHM <- function () 
{
  return(format(Sys.time(), "%H:%M"))
}
