#' @title viral_lifestyle
#' @description This function assigns the lifestyle, lytic or temperate/lysogenic, to the virus with known host
#' #'
#' @param vcontact_df dataframe from vcontact2
#' @import dplyr
#' @import janitor
#' @examples
#' \dontrun{
#' quiet(cat("test"))
#' quiet(warning("test"))
#' quiet(warning("test"), all=T)
#' }
#' # This is a function that suppresses log info
#' @export
viral_lifestyle <- function(vcontact_df) {
  '%>%' <- tidyr::`%>%`
  v_lifestyle<-readRDS(system.file("External_files", "v_lifestyle.rds", package = "MetaViral", mustWork = TRUE))

  vcontact_df$Lifestyle <- NULL

  lifestyle<-sapply(seq_along(vcontact_df$Genome), function(x){

    if(janitor::make_clean_names(vcontact_df$Genome[x]) %in% v_lifestyle$virus_name_c==TRUE){
      temp<-which(v_lifestyle$virus_name_c==janitor::make_clean_names(vcontact_df$Genome[x])[1])
      temp2<-v_lifestyle$lifestyle[temp]

    }else if(janitor::make_clean_names(vcontact_df$Genome[x]) %in% v_lifestyle$host_name_c==TRUE){
      temp<-which(v_lifestyle$host_name_c==janitor::make_clean_names(vcontact_df$Genome[x])[1])
      temp2<-v_lifestyle$lifestyle[temp]

    }else{
      temp2<-NA
    }

    return(temp2)
  })#end of sapply

    vcontact_df2<-cbind(lifestyle, vcontact_df)

  return(vcontact_df2)

}#end of function host_domain

