#' @title cleaning
#' @description cleans
#' #'
#' @param data
#' @param sepby regex comands to separate the column. Vcontact2 regex: "-(?=[^-])" ; DRAM-v regex: "_(?=[^_]+$)"
#' @param fullclean
#' @param output_from
#' @param ...
#' @import dplyr
#' @importFrom tidyr %>%
#' @examples
#' \dontrun{
#' quiet(cat("test"))
#' quiet(warning("test"))
#' quiet(warning("test"), all=T)
#' }
#' # This is a function that suppresses log info
#' @export
cleaning<-function(data, sepby="-(?=[^-])", fullclean=TRUE, output_from="vcontact2", ...){
  '%>%' <- tidyr::`%>%`
  if (fullclean==TRUE){
    viral_genome <- tidyr::separate(data, col=filename, into=c("Biome", "Category"), sep=sepby, extra = "merge")
    viral_genome <- tidyr::separate(viral_genome, col=Biome, into=c("junk", "Biome"), sep="(?=[A-z0-9])", extra = "merge")
    viral_genome <- tidyr::separate(viral_genome, col=Category, into=c("junk", "Category"), sep="(?=[0-9])", remove = FALSE)
    viral_genome <- tidyr::separate(viral_genome, col=Category, into=c("Category", "junk"), sep="[-]", extra = "merge")
    viral_genome$junk <- NULL
    viral_genome<-viral_genome[,-1]

  }else {
    viral_genome <- tidyr::separate(data, col=filename, into=c("Biome", "Category"), sep="-(?=[^-])", extra = "merge")
  }

  if (output_from=="vcontact2"){
    #sort and group by each VC
    viral_vc<-viral_genome %>%
      dplyr::group_by(`VC Subcluster`) %>%
      dplyr::arrange(Genome, Family, Order, Genus, Biome)

    viral_vc<-viral_vc[!is.na(viral_vc$VC), ]

    viral_vc <- tidyr::separate(viral_vc, col=VC, into=c("VC", "VC2"), sep="\\D", extra = "merge")

    viral_vc<-transform(viral_vc, VC = as.numeric(VC))
    viral_vc<-transform(viral_vc, VC2 = as.numeric(VC2))

    #viral_sorted <- viral_vc[order(as.integer(viral_vc$VC),decreasing = FALSE), ]

    viral_vc <- viral_vc[order(viral_vc$VC, viral_vc$VC2, decreasing = FALSE), ]

    return(viral_vc)

  }else if (output_from=="DRAMv") {

    return(viral_annotations)

  }else {
    warning("output_from must be vcontact2 or DRAMv")
  }
}



