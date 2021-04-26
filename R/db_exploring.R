#' @title db_exploring
#' @description isolates the database selected, counts for the biome and can perform other manipulations
#' #'
#' @param viral_annotations
#' @param database
#' @import dplyr
#' @import KEGGREST
#' @examples
#' \dontrun{
#' quiet(cat("test"))
#' quiet(warning("test"))
#' quiet(warning("test"), all=T)
#' }
#' # This is a function that suppresses log info
#' @export
db_exploring<-function(viral_annotations, database="all"){
  '%>%' <- tidyr::`%>%`
  if (database=="peptidase"){
    peptidase<-dplyr::select(viral_annotations, peptidase_hit, peptidase_id, peptidase_family, peptidase_RBH, peptidase_eVal, Biome, fasta)
    peptidase<-na.omit(peptidase)
    peptidase<-peptidase %>%
      dplyr::add_count(peptidase_hit,name="peptidase_counts", Biome)
    return(peptidase)

  }else if (database=="vog"){
    vog<-dplyr::select(viral_annotations, vogdb_description, Biome, fasta)
    vog<-na.omit(vog)
    vog<-vog %>%
      dplyr::add_count(vogdb_description,name="vog_counts", Biome)
    return(vog)

  }else if (database=="kegg"){
    kegg<-dplyr::select(viral_annotations, Biome, kegg_hit, kegg_id, fasta)
    kegg<-na.omit(kegg)
    kegg<-kegg %>%
      dplyr::add_count(kegg_hit,name="kegg_counts", Biome)

    kegg$Pathway<-NA
    kegg$Pathway2<-NA
    kegg$Pathway3<-NA

    for (i in 1:nrow(kegg)) {
      protein_id<-c(kegg$kegg_id[i])
      query<-tryCatch(KEGGREST::keggGet(protein_id), error=function(e) NULL)

      if (length(query[[1]]$PATHWAY)==1) {
        kegg$Pathway[i]<-query[[1]]$PATHWAY

      } else if (length(query[[1]]$PATHWAY)==2) {
        kegg$Pathway[i]<-query[[1]]$PATHWAY[1]
        kegg$Pathway2[i]<-query[[1]]$PATHWAY[2]

      } else if (length(query[[1]]$PATHWAY)==3) {
        kegg$Pathway[i]<-query[[1]]$PATHWAY[1]
        kegg$Pathway2[i]<-query[[1]]$PATHWAY[2]
        kegg$Pathway3[i]<-query[[1]]$PATHWAY[3]
      }
    }#end of for loop

    kegg<-kegg %>%
      dplyr::add_count(Pathway,name="pathway_counts", Biome)

    return(kegg)

  }else if (database=="pfam"){
    pfam<-dplyr::select(viral_annotations, Biome, pfam_hits, pfam_id, fasta)
  pfam<-na.omit(pfam)
  pfam<-pfam %>%
    dplyr::add_count(pfam_hits,name="pfam_counts", Biome)
  return(pfam)

  }else if (database=="viral"){
    viraldb<-dplyr::select(viral_annotations, Biome, viral_function, viral_id, viral_tax, viral_identity,
    viral_bitScore, viral_eVal, fasta)
  viraldb<-na.omit(viraldb)
  viraldb<-viraldb %>%
    dplyr::add_count(viral_function,name="viral_counts", Biome)
  return(viraldb)

  }else if (database=="all"){
    peptidase<-dplyr::select(viral_annotations, peptidase_hit, peptidase_id, peptidase_family, peptidase_RBH, peptidase_eVal, Biome)
    peptidase<-na.omit(peptidase)
    peptidase<-peptidase %>%
      dplyr::add_count(peptidase_hit,name="peptidase_counts", Biome)

    vog<-dplyr::select(viral_annotations, vogdb_description, Biome)
    vog<-na.omit(vog)
    vog<-vog %>%
      dplyr::add_count(vogdb_description,name="vog_counts", Biome)

    kegg<-dplyr::select(viral_annotations, Biome, kegg_hit, kegg_id)
    kegg<-na.omit(kegg)
    kegg<-kegg %>%
      dplyr::add_count(kegg_hit,name="kegg_counts", Biome)

    kegg$Pathway<-NA
    kegg$Pathway2<-NA
    kegg$Pathway3<-NA

    for (i in 1:nrow(kegg)) {
      protein_id<-c(kegg$kegg_id[i])
      query<-tryCatch(KEGGREST::keggGet(protein_id), error=function(e) NULL)

      if (length(query[[1]]$PATHWAY)==1) {
        kegg$Pathway[i]<-query[[1]]$PATHWAY

      } else if (length(query[[1]]$PATHWAY)==2) {
        kegg$Pathway[i]<-query[[1]]$PATHWAY[1]
        kegg$Pathway2[i]<-query[[1]]$PATHWAY[2]

      } else if (length(query[[1]]$PATHWAY)==3) {
        kegg$Pathway[i]<-query[[1]]$PATHWAY[1]
        kegg$Pathway2[i]<-query[[1]]$PATHWAY[2]
        kegg$Pathway3[i]<-query[[1]]$PATHWAY[3]
      }
    }#end of for loop

    kegg<-kegg %>%
      dplyr::add_count(Pathway,name="pathway_counts", Biome)

    pfam<-dplyr::select(viral_annotations, Biome, pfam_hits, pfam_id)
    pfam<-na.omit(pfam)
    pfam<-pfam %>%
      dplyr::add_count(pfam_hits,name="pfam_counts", Biome)

    viraldb<-dplyr::select(viral_annotations, Biome, viral_function, viral_id, viral_tax, viral_identity,
                           viral_bitScore, viral_eVal)
    viraldb<-na.omit(viraldb)
    viraldb<-viraldb %>%
      dplyr::add_count(viral_function,name="viral_counts", Biome)
    #merge all db
    database_all<-Reduce(function(x,y) merge(x,y,by="Biome",all=TRUE) ,list(kegg, pfam, vog, viraldb, peptidase))
    return(database_all)

  }else {
    warning("choose one of the following: kegg, pfam, vog, viraldb, peptidase or all")
  }
}#end of function
