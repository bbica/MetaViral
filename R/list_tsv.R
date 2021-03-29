#' @title list_tsv
#' @description This is a function that suppresses log info
#' #'
#' @param flnm name of the file
#' @importFrom readr read_tsv
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
list_tsv <- function(flnm) {
  '%>%' <- tidyr::`%>%`
  return_namefile <-
    readr::read_tsv(flnm) %>%
    dplyr::mutate(filename = flnm)
  return(return_namefile)
}