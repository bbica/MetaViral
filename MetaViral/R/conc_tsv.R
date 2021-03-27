#' @title conc_tsv
#' @description This is a function that suppresses log info
#' #'
#' @param x the function to suppress the log info
#' @param all If TRUE then suppress warnings and messages as well; otherwise, only suppress printed output (such as from print or cat).
#' @importFrom tidyr %>%
#' @importFrom purrr map_df
#' @examples
#' \dontrun{
#' quiet(cat("test"))
#' quiet(warning("test"))
#' quiet(warning("test"), all=T)
#' }
#' # This is a function that suppresses log info
#' @export

conc_tsv <- function(path){
  '%>%' <- tidyr::`%>%`
  wd_dram<-path
  setwd(wd_dram)
  export_tsv <-
    list.files(pattern = "*.tsv",
               full.names = T) %>%
    purrr::map_df(~list_tsv(.))
  return(export_tsv)
}


