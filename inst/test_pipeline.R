library(MetaViral)

# Clear R's brain
rm(list = ls())


#vcontact2
wd_genomes<-"~/Bernardo/ALL DATA 2020/Vcontact2/genome_by_genome/"
setwd(wd_genomes)

viral_genome<-import_files("~/Bernardo/ALL DATA 2020/Vcontact2/genome_by_genome/", filetype = "csv")

#Cleaning...
viral_genome <- tidyr::separate(viral_genome, col=filename, into=c("Biome", "Category"), sep="-(?=[^-])", extra = "merge")
viral_genome <- tidyr::separate(viral_genome, col=Biome, into=c("junk", "Biome"), sep="(?=[A-z0-9])", extra = "merge")
viral_genome <- tidyr::separate(viral_genome, col=Category, into=c("junk", "Category"), sep="(?=[0-9])", remove = FALSE)
viral_genome <- tidyr::separate(viral_genome, col=Category, into=c("Category", "junk"), sep="[-]", extra = "merge")
viral_genome$junk <- NULL
viral_genome<-viral_genome[,-1]


viral_genome_mx<-abundance_standart(viral_genome)

#VC

#sort and group by each VC
viral_vc<-viral_genome %>%
  group_by(`VC Subcluster`) %>%
  arrange(Genome, Family, Order, Genus, Biome, Category)

viral_vc<-viral_vc[!is.na(viral_vc$VC), ]

viral_vc <- tidyr::separate(viral_vc, col=VC, into=c("VC", "VC2"), sep="\\D", extra = "merge")

viral_vc<-transform(viral_vc, VC = as.numeric(VC))
viral_vc<-transform(viral_vc, VC2 = as.numeric(VC2))

#viral_sorted <- viral_vc[order(as.integer(viral_vc$VC),decreasing = FALSE), ]

viral_vc <- viral_vc[order(viral_vc$VC, viral_vc$VC2, decreasing = FALSE), ]

#abundance_vc function
falta_mx<-MetaViral::abundance_vc(viral_vc=viral_vc, abuntype = "relative", taxa = "Genome",output_type = "matrix" )



#DRAMv

viral_genome<-import_files("~/Bernardo/ALL DATA 2020/DRAMv/annotations/", filetype = "tsv")

