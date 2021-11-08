#' @title abundance_standard
#' @description This is a function that suppresses log info
#' #'
#' @param viral_df dataframe, gene_to_genome.csv, output from VContact2
#' @param taxa taxon selected by the user among the following: "Genome","Order","Family","Genus". Default set to "Family"
#' @param abuntype abundance to construct the matrix. "relative" or "absolute"; default set to "relative"
#' @importFrom purrr as_vector
#' @examples
#' \dontrun{
#' quiet(cat("test"))
#' quiet(warning("test"))
#' quiet(warning("test"), all=T)
#' }
#' # This is a function that suppresses log info
#' @export
abundance_standard<-function(viral_df, taxa="Family", abuntype="relative"){
  biomes<-unique(viral_df$Biome)
  number_biomes<-length(biomes)

  #if statement for several taxas
  vector_columns<-colnames(viral_df)
  vector_taxas<-c("Genome", "Order","Family","Genus")
  correspondent_taxas<-which(vector_columns%in%vector_taxas)

  if(taxa==vector_taxas[3]){ #in this case, Family
    column_number<-which(vector_columns==vector_taxas[3])
    taxa_list<-unique(viral_df[,column_number])
    number_taxa<-length(taxa_list)

  } else if (taxa==vector_taxas[1]){
    column_number<-which(vector_columns==vector_taxas[1])# in case of Genome
    taxa_list<-unique(viral_df[,column_number])
    number_taxa<-length(taxa_list)

  } else if (taxa==vector_taxas[2]){
    column_number<-which(vector_columns==vector_taxas[2])# in case of Order
    taxa_list<-unique(viral_df[,column_number])
    number_taxa<-length(taxa_list)

  } else if (taxa==vector_taxas[4]){
    column_number<-which(vector_columns==vector_taxas[4])# in case of Genus
    taxa_list<-unique(viral_df[,column_number])
    number_taxa<-length(taxa_list)

  } else {
    warning(paste0("Taxa selected isn't among the following: ", paste(vector_taxas, collapse=', ')))
  }

  #Construction  of a "hollow" matrix (taxa_matrix)with the biomes as colum names and the taxas as rownames.

  taxa_matrix<-matrix(data=NA, nrow=number_taxa, ncol=number_biomes)
  colnames(taxa_matrix)<-biomes
  rownames(taxa_matrix)<-purrr::as_vector(taxa_list)


  for (i in 1:length(biomes)){
    bioma<-biomes[i]
    temp<-which(viral_df$Biome==bioma)
    subset<-viral_df[temp,]
    dat<-subset[,column_number]

    if (abuntype=="relative"){
      abun_type<-as.matrix(table(dat)/nrow(table(dat)))

    } else if (abuntype=="absolute"){
      abun_type<-as.matrix(table(dat))
    }

    subset_col_viral_df<-cbind(rownames(abun_type),abun_type[,1])
    rownames(subset_col_viral_df)<-NULL
    subset_col_viral_df_order<-subset_col_viral_df[order(match(subset_col_viral_df[,1], abun_type)),]
    colnames(subset_col_viral_df_order)[2]<-"relAbundByBiome"

    for (j in 1:length(purrr::as_vector(taxa_list))){
      temp_taxa<-purrr::as_vector(taxa_list[j])
      taxa_selected<-which(subset_col_viral_df_order[,1]==temp_taxa)
      abun_of_taxa<-subset_col_viral_df_order[taxa_selected,2] #gives abundance of the taxa selected
      abun_of_taxa<-as.numeric(abun_of_taxa)
      #temp2<-which(subset[,column_number]==temp_taxa)
      #mini_subset<-subset[temp2,]
      #soma_abun<-sum(mini_subset$relAbundByBiome)
      taxa_matrix[j,i]<-abun_of_taxa
    }
  }
  return(taxa_matrix)
}
