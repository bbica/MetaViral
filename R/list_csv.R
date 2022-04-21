#' @title list_csv
#' @description This is a function that suppresses log info
#' #'
#' @param flnm name of the file
#' @param col_names3 header option, default set to TRUE
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
list_csv <- function(flnm, col_names3=TRUE) {
  '%>%' <- tidyr::`%>%`
  if (col_names3==TRUE){
    return_namefile <-
      data.table::fread(flnm, header = TRUE, sep=",") %>%
      dplyr::mutate(filename = flnm)
    return(return_namefile)

  }else{
    return_namefile <-
      data.table::fread(flnm, header =FALSE, sep=",") %>%
      dplyr::mutate(filename = flnm)
    return(return_namefile)
  }
}

