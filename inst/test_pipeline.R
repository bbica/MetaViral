library(MetaViral)

# Clear R's brain
rm(list = ls())

#VCONTACT2
#set working directory to where the genome_by_genome_overview.csv files are located
wd_dir_exe<-"C:/Users/bbica/Desktop/MetaViral/inst/Example_Data"
setwd(wd_dir_exe)

wd_genomes<-"~/Bernardo/ALL DATA 2020/Vcontact2/genome_by_genome/"
setwd(wd_genomes)

#Import with import_files function that calls the list_csv function (filetype="csv")
viral_genome<-import_files(wd_genomes, filetype = "csv")

#Cleaning...
viral_vc<-cleaning(viral_genome, output_from="vcontact2")

#Creation of a abundance matrix with abundance_vc function
VC_mx_abs<-MetaViral::abundance_vc(viral_vc=viral_vc, abuntype = "absolute", taxa = "Genome",output_type = "matrix" )
VC_mx_rel<-MetaViral::abundance_vc(viral_vc=viral_vc, abuntype = "relative", taxa = "Genome",output_type = "matrix" )

#bio_statistics function
VC_stats<-bio_statistics(VC_mx_abs)


#DRAMv

viral_genome<-import_files("~/Bernardo/ALL DATA 2020/DRAMv/annotations/", filetype = "tsv")

