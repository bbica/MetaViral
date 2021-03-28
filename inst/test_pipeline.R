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

#DRAMv

viral_genome<-import_files("~/Bernardo/ALL DATA 2020/DRAMv/annotations/", filetype = "tsv")




