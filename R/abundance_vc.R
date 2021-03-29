#' @title abundance_vc
#' @description This is a function that suppresses log info
#' @param viral_vc dataframe, gene_to_genome.csv, output from VContact2
#' @param taxa taxon selected by the user among the following: "Genome","Order","Family","Genus". Default set to "Family"
#' @param abuntype abundance to construct the matrix. "relative" or "absolute"; default set to "relative"
#' @importFrom purrr as_vector
#' @importFrom dplyr add_count group_by mutate select
#' @examples
#' \dontrun{
#' quiet(cat("test"))
#' quiet(warning("test"))
#' quiet(warning("test"), all=T)
#' }
#' # This is a function that suppresses log info
#' @export
abundance_vc<-function(viral_vc, taxa="Family", abuntype="relative"){
#Majority rules
viral_test1<-viral_vc%>%
  # add a column n with count by categories
  add_count(VC.Subcluster, Genus, Biome, Category) %>%
  # select max or first occurrence
  group_by(VC.Subcluster) %>%
  # keep only first TRUE
  mutate(Majority_g = Genus[n == max(n)][1]) %>%
  # do not keep temp var
  select(-n)

viral_test1<-viral_test1 %>%
  # add a column n with count by categories
  add_count(VC.Subcluster, Family, Biome, Category) %>%
  # select max or first occurrence
  group_by(VC.Subcluster) %>%
  # keep only first TRUE
  mutate(Majority_f = Family[n == max(n)][1]) %>%
  # do not keep temp var
  select(-n)

viral_test1<-viral_test1 %>%
  # add a column n with count by categories
  add_count(VC.Subcluster, Genome, Biome, Category) %>%
  # select max or first occurrence
  group_by(VC.Subcluster) %>%
  # keep only first TRUE
  mutate(Majority_sp = Genome[n == max(n)][1]) %>%
  # do not keep temp var
  select(-n)

viral_test1<-viral_test1 %>%
  # add a column n with count by categories
  add_count(VC.Subcluster, Order, Biome, Category) %>%
  # select max or first occurrence
  group_by(VC.Subcluster) %>%
  # keep only first TRUE
  mutate(Majority_o = Order[n == max(n)][1]) %>%
  # do not keep temp var
  select(-n)


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
          Majority_taxa<-"Majority_f"
        }else if (taxa=="Order"){
          Majority_taxa<-"Majority_o"
        }else if (taxa=="Genus"){
          Majority_taxa<-"Majority_g"
        }else if (taxa=="Genome"){
          Majority_taxa<-"Majority_sp"
        }else {
          Majority_taxa<-"Majority_f"
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
return(vec4)
}#end of function