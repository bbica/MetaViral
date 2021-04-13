#' @title abundance_vc
#' @description Creates a abundance matrix from the genome_by_genome_overview.csv files outputed by VCONTACT2
#' @param viral_vc dataframe, gene_to_genome.csv, output from VContact2
#' @param taxa taxon selected by the user among the following: "Genome","Order","Family","Genus". Default set to "Family"
#' @param abuntype abundance to construct the matrix. "relative" or "absolute"; default set to "relative"
#' @importFrom purrr as_vector
#' @import dplyr
#' @importFrom rlang is_empty
#' @examples
#' \dontrun{
#' quiet(cat("test"))
#' quiet(warning("test"))
#' quiet(warning("test"), all=T)
#' }
#' # This is a function that suppresses log info
#' @export
abundance_vc<-function(viral_vc, taxa="Family",output_type="matrix", abuntype="relative"){
  '%>%' <- tidyr::`%>%`
  #Majority rules
  viral_test1<-viral_vc %>%
    # add a column n with count by categories
    dplyr::add_count(VC.Subcluster, Genus, Biome, Category) %>%
    # select max or first occurrence
    dplyr::group_by(VC.Subcluster) %>%
    # keep only first TRUE
    dplyr::mutate(VC_Genus = Genus[n == max(n)][1]) %>%
    # do not keep temp var
    dplyr::select(-n)

  viral_test1<-viral_test1 %>%
    dplyr::add_count(VC.Subcluster, Family, Biome, Category) %>%
    dplyr::group_by(VC.Subcluster) %>%
    dplyr::mutate(VC_Family = Family[n == max(n)][1]) %>%
    dplyr::select(-n)

  viral_test1<-viral_test1 %>%
    dplyr::add_count(VC.Subcluster, Genome, Biome, Category) %>%
    dplyr::group_by(VC.Subcluster) %>%
    dplyr::mutate(VC_Genome = Genome[n == max(n)][1]) %>%
    dplyr::select(-n)

  viral_test1<-viral_test1 %>%
    dplyr::add_count(VC.Subcluster, Order, Biome, Category) %>%
    dplyr::group_by(VC.Subcluster) %>%
    dplyr::mutate(VC_Order = Order[n == max(n)][1]) %>%
    dplyr::select(-n)


  biomas<-unique(viral_vc$Biome)
  categories<-unique(viral_vc$Category)
  viralclusters<-unique(viral_vc$VC.Subcluster)
  vec<-NULL
  for (i in 1:length(biomas)){
    for (j in 1:length(categories)){
      for (k in 1:length(viralclusters)){
        rowselected<-which(viral_test1$Biome==biomas[i]&viral_test1$Category==categories[j]&viral_test1$VC.Subcluster==viralclusters[k])
        vector_id<-paste0(biomas[i], "-", categories[j], "-", viralclusters[k])
        occurencies<-length(rowselected)
        if (taxa=="Family"){
          Majority_taxa<-"VC_Family"
        }else if (taxa=="Order"){
          Majority_taxa<-"VC_Order"
        }else if (taxa=="Genus"){
          Majority_taxa<-"VC_Genus"
        }else if (taxa=="Genome"){
          Majority_taxa<-"VC_Genome"
        }else {
          Majority_taxa<-"VC_Family"
          print("Family selected as the taxa")
        }

        taxa_column<-which(colnames(viral_test1)==Majority_taxa)
        correspondent_taxa<-viral_test1[rowselected[1], taxa_column]
        vec<-rbind(vec, cbind(vector_id, correspondent_taxa, occurencies))
      }#loop k
    }#loop j
    print(biomas[i])
  }#loop i
  vec2<-vec[!is.na(vec[,2]), ]
  vec_biomas<-c()
  for (i in 1:nrow(vec2)){
    vec_biomas<-c(vec_biomas,strsplit(vec2[i,1], split = "-")[[1]][1])
  }
  vec_biomas_mx<-c()
  for (i in 1:length(biomas)){
    assign(paste0(biomas[i],"_data"),vec2[which(vec_biomas==biomas[i]),])
    vec_biomas_mx<-c(vec_biomas_mx, paste0(biomas[i],"_data"))
  }
  vec4<-list()
  for (i in 1:length(vec_biomas_mx)){
    bioma_mx<-get(vec_biomas_mx[i])
    unique_taxa<-unique(bioma_mx[,2])
    vec3<-NULL
    for (j in 1:length(unique_taxa)){
      if (abuntype=="absolute"){
        sum_abs<-sum(bioma_mx[which(bioma_mx[,2]==unique_taxa[j]), 3])
        vec3<-rbind(vec3,cbind(unique_taxa[j], sum_abs))

      }else{
        sum_abs<-sum(bioma_mx[which(bioma_mx[,2]==unique_taxa[j]), 3])
        sum_rel<-sum_abs/sum(bioma_mx[,3])
        vec3<-rbind(vec3,cbind(unique_taxa[j], sum_rel))
      }
      #dividor<-nrow(vec3)/length(vec_biomas_mx)
      #vec4<-vec3[1:dividor,]
    }
    vec4[[i]]<-vec3
  }
  names(vec4)<-biomas
  if (output_type=="list"){
    return(vec4)

  }else{
    vecall<-c()

    for (i in 1:length(vec4)){
      vectmp<-vec4[[i]][,1]
      vecall<-c(vecall,vectmp)
    }
    vecuniquetaxa<-unique(vecall)
    matrix_final<-matrix(data=NA, nrow=length(biomas), ncol=length(vecuniquetaxa))
    rownames(matrix_final)<-biomas
    colnames(matrix_final)<-vecuniquetaxa
    for (i in 1:length(biomas)){
      for (j in 1:length(vecuniquetaxa)){
        corresp<-which(vec4[[i]][,1]==vecuniquetaxa[j])
        insert_in_matrix<-c()
        if (rlang::is_empty(corresp)==TRUE){
          insert_in_matrix<-0
        }else{
          insert_in_matrix<-as.numeric(vec4[[i]][corresp,2])
        }
        matrix_final[i,j]<-insert_in_matrix
      }
    }
    return(matrix_final)
  }#end else
}#end of function
