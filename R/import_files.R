#' @title import_files
#' @description imports files
#' #'
#' @param path path where the files are present
#' @param filetype example: "tsv" "csv". Default set to tsv
#' @import dplyr
#' @importFrom purrr map_df
#' @examples
#' \dontrun{
#' quiet(cat("test"))
#' quiet(warning("test"))
#' quiet(warning("test"), all=T)
#' }
#' # This is a function that suppresses log info
#' @export

import_files <- function(path, filetype="tsv"){
  '%>%' <- tidyr::`%>%`
  wd_dram<-path
  setwd(wd_dram)
  if (filetype=="tsv"){
    export_files <-
      list.files(pattern = paste0("*.", filetype),
                 full.names = T) %>%
      purrr::map_df(~list_tsv(.))
  } else if (filetype=="csv"){
    export_files <-
      list.files(pattern = paste0("*.", filetype),
                 full.names = T) %>%
      purrr::map_df(~list_csv(.))
  } else {
    warning("Filetype must be 'csv' or 'tsv'")
    }

  return(export_files)
}

#end
