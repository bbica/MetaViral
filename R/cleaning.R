#' @title cleaning
#' @description cleans the output from vcontact2 or DRAM-v
#' #'
#' @param viral_data
#' @param output_from
#' @param ...
#' @import dplyr
#' @importFrom gtools permutations
#' @examples
#' \dontrun{
#' quiet(cat("test"))
#' quiet(warning("test"))
#' quiet(warning("test"), all=T)
#' }
#' # This is a function that suppresses log info
#' @export
cleaning<-function(viral_data, output_from="vcontact2", ...){
  '%>%' <- tidyr::`%>%`
  if (output_from=="vcontact2"){
    viral_genome <- tidyr::separate(viral_data, col=filename, into=c("Biome", "Category"), sep="-(?=[^-])", extra = "merge")
    viral_genome <- tidyr::separate(viral_genome, col=Biome, into=c("junk", "Biome"), sep="(?=[A-z0-9])", extra = "merge")
    viral_genome <- tidyr::separate(viral_genome, col=Category, into=c("junk", "Category"), sep="(?=[0-9])", remove = FALSE)
    viral_genome <- tidyr::separate(viral_genome, col=Category, into=c("Category", "junk"), sep="[-]", extra = "merge")
    viral_genome$junk <- NULL
    viral_genome<-viral_genome[,-1]

    viral_vc<-viral_genome %>%
      dplyr::group_by(`VC Subcluster`) %>%
      dplyr::arrange(Genome, Family, Order, Genus, Biome)

    viral_vc<-viral_vc[!is.na(viral_vc$VC), ]

    viral_vc <- tidyr::separate(viral_vc, col=VC, into=c("VC", "VC2"), sep="\\D", extra = "merge")

    viral_vc<-transform(viral_vc, VC = as.numeric(VC))
    viral_vc<-transform(viral_vc, VC2 = as.numeric(VC2))
    viral_vc <- viral_vc[order(viral_vc$VC, viral_vc$VC2, decreasing = FALSE), ]

    return(viral_vc)



  }else if (output_from=="DRAMv") {
    #extract only the name of the biome from the filename column (FAZER pipes!!!)
    DRAMv_with_sources <- tidyr::separate(viral_data, col=filename, into=c("Biome", "junk"), sep="_(?=[^_]+$)", extra = "merge")
    DRAMv_with_sources <- tidyr::separate(DRAMv_with_sources, col=Biome, into=c("junk2", "Biome"), sep="(?=[A-z0-9])", extra = "merge")
    DRAMv_with_sources$junk <- NULL
    DRAMv_with_sources$junk2 <- NULL

    DRAMv_with_sources2 <- DRAMv_with_sources%>%
      tidyr::separate(col=viral_hit, into=c("viral_function", "protein_and_vname"), sep="\\s", extra="merge")%>%
      tidyr::separate(col=protein_and_vname, into=c("viral_function", "viral_tax"), sep="(?=\\[)")%>% #separates the column by the fist  "["
      tidyr::separate(col=kegg_hit, into=c("kegg_hit", "kegg_hit2", "kegg_hit3", "kegg_hit4"), sep=";")%>%
      tidyr::separate(col=kegg_hit, into=c("kegg_hit", "junk"), sep="(?=\\[)")%>%
      tidyr::separate(col=vogdb_description, into=c("vogdb_description", "junk"), sep=";")%>%
      tidyr::separate(col=vogdb_description, into=c("junk2", "vogdb_description"), sep="(?=\\s)", extra = "merge")%>%
      tidyr::separate(col=pfam_hits, into=c("pfam_hits", "pfam_id"), sep="(?=\\[)")

    DRAMv_with_sources2$junk <- NULL
    DRAMv_with_sources2$junk2 <- NULL

    #Extracts the taxonomic predition from the square brackets
    DRAMv_with_sources2$viral_tax = sub(".*\\[([^][]+)].*", "\\1", DRAMv_with_sources2$viral_tax)

    #removes rows without any predictions
    DRAMv_with_sources3 <- DRAMv_with_sources2 %>%
      dplyr::filter_at( dplyr::vars(kegg_hit,viral_function,pfam_hits,cazy_hits,vogdb_description,peptidase_hit),any_vars(!is.na(.)))

    #removes rows with NA in the amg_flags
    DRAMv_with_sources4 <- DRAMv_with_sources3[!is.na(DRAMv_with_sources3$amg_flags),]
    #filter for only auxiliary scores 4 and 5
    DRAMv_with_sources4<-dplyr::filter(DRAMv_with_sources4, auxiliary_score<4)
    #filter out flags related with viral function
    unwanted_flags <- c("F", "T", "P", "A")

      for (i in 1:max(nchar(DRAMv_with_sources4$amg_flags))) {
        prm<-gtools::permutations(n=4,r=i,v=unwanted_flags, repeats.allowed = TRUE)
        prm<-apply(prm, 1, function(x) paste0(x, collapse=''))
        DRAMv_with_sources4<-DRAMv_with_sources4[ ! DRAMv_with_sources4$amg_flags %in% prm, ]
      }

    viral_annotations<-DRAMv_with_sources4
    return(viral_annotations)

  }else {
    warning("output_from must be vcontact2 or DRAMv")
  }
}



