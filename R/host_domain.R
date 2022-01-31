#' @title host_domain
#' @description This function assigns to the host, or genome, the respective domain, archeae, bacteria or eukarya
#' #'
#' @param vcontact_df dataframe from vcontact2
#' @param reduce_names if TRUE separates the Genome column so that only the host genus is present, default set to FALSE
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
host_domain <- function(vcontact_df, reduce_names=FALSE) {
  '%>%' <- tidyr::`%>%`
  hostdb<-readRDS(system.file("External_files", "hostdb.rds", package = "MetaViral", mustWork = TRUE))
  vcontact_df2 <- tidyr::separate(vcontact_df, col=Genome, into=c("Genome", "Genome2"), sep="[[:punct:]]", extra = "merge")
  vcontact_df2$Genome2 <- NULL


  domain<-sapply(seq_along(vcontact_df$Genome), function(x){

    if(janitor::make_clean_names(vcontact_df2$Genome[x]) %in% hostdb$Genus2==TRUE){
      temp<-which(hostdb$Genus2==janitor::make_clean_names(vcontact_df2$Genome[x])[1])
      temp2<-hostdb$Domain[temp]

    }else if(janitor::make_clean_names(vcontact_df2$Genome[x]) %in% hostdb$Genus==TRUE){
      temp<-which(hostdb$Genus==janitor::make_clean_names(vcontact_df2$Genome[x])[1])
      temp2<-hostdb$Domain[temp]

    }else if(janitor::make_clean_names(vcontact_df2$Genome[x]) %in% hostdb$virus_name==TRUE){
      temp<-which(hostdb$virus_name==janitor::make_clean_names(vcontact_df2$Genome[x])[1])
      temp2<-hostdb$Domain[temp]

    }else{
      temp2<-NA
    }

    return(temp2)
  })#end of sapply

  if(reduce_names==TRUE){
    vcontact_df3<-cbind(domain, vcontact_df2)
    vcontact_df3$Genome <- paste0(vcontact_df3$Genome, "-",  vcontact_df3$domain)

  }else{
    vcontact_df3<-cbind(domain, vcontact_df)
    vcontact_df3$Genome <- paste0(vcontact_df3$Genome, "-",  vcontact_df3$domain)
  }

  return(vcontact_df3)

}#end of function host_domain

