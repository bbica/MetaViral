#' @title bio_statistics
#' @description performs multiple diversity indexes
#' #'
#' @param abundance_mx
#' @import vegan
#' @examples
#' \dontrun{
#' quiet(cat("test"))
#' quiet(warning("test"))
#' quiet(warning("test"), all=T)
#' }
#' # This is a function that suppresses log info
#' @export
bio_statistics<-function(abundance_mx){
  Biomes<-rownames(abundance_mx)
  df_stats <- as.data.frame(matrix(0, ncol = 1, nrow = length(Biomes)))
  df_stats[,1]<-Biomes
  colnames(df_stats)<-"Biomes"

  df_stats$Abundance<-apply(abundance_mx,1,sum)#Total number of individuals(Abundance)
  df_stats$Richness<-apply(abundance_mx>0,1,sum)#Taxa Richness
  df_stats$Menhinick<-df_stats$Richness/sqrt(df_stats$Abundance)#Menhinick's index
  df_stats$Margalef<-(df_stats$Richness-1)/log(df_stats$Abundance)#Margalef's index
  df_stats$Shannon<-vegan::diversity(abundance_mx, index="shannon")#Shannon-Wiener Index
  df_stats$Simpson<-vegan::diversity(abundance_mx, index="simpson")#Simpson's Index

  df_stats$Pilou_evenness<-df_stats$Simpson/log(df_stats$Richness)#Pilou evenness
  df_stats$Simpson_evenness<-exp(df_stats$Simpson)/df_stats$Richness#Hill's evenness ratio for Simpson's Index
  df_stats$Shannon_evenness<-exp(df_stats$Shannon)/df_stats$Richness#Hill's evenness ratio for Shannon-Wiener Index

  df_stats$True_Shannon<-exp(df_stats$Shannon)#True diversity for Shannon-Wiener Index
  df_stats$True_Simpson<-1/(df_stats$Simpson)#True diversity for Simpson's Index

  return(df_stats)
}#end of function





