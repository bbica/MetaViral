#' @title list_tsv
#' @description This is a function that suppresses log info
#' #'
#' @param flnm name of the file
#' @param col_names2 TRUE if first row is the header, default
#' @importFrom data.table fread
#' @import dplyr
#' @examples
#' \dontrun{
#' quiet(cat("test"))
#' quiet(warning("test"))
#' quiet(warning("test"), all=T)
#' }
#' # This is a function that suppresses log info
#' @export
list_tsv <- function(flnm, col_names2=TRUE) {
  '%>%' <- tidyr::`%>%`

  if (col_names2==TRUE){
    return_namefile <-
      data.table::fread(flnm, header = TRUE, sep="\t") %>%
      dplyr::mutate(filename = flnm)
    return(return_namefile)

  }else{
    return_namefile <-
      data.table::fread(flnm, header = FALSE, sep="\t") %>%
      dplyr::mutate(filename = flnm)
    return(return_namefile)

  }
}
