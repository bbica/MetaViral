#' @title list_csv
#' @description This is a function that suppresses log info
#' #'
#' @param flnm name of the file
#' @importFrom readr read_csv
#' @importFrom dplyr mutate
#' @importFrom tidyr %>%
#' @examples
#' \dontrun{
#' quiet(cat("test"))
#' quiet(warning("test"))
#' quiet(warning("test"), all=T)
#' }
#' # This is a function that suppresses log info
#' @export
list_csv <- function(flnm) {
  '%>%' <- tidyr::`%>%`
  return_namefile <-
    readr::read_csv(flnm) %>%
    dplyr::mutate(filename = flnm)
  return(return_namefile)
}
