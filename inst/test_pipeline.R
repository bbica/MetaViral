library(MetaViral)

# Clear R's brain
rm(list = ls())

#VCONTACT2
#set working directory to where the genome_by_genome_overview.csv files are located
wd_example<-system.file("Example_Data", package= "MetaViral")
setwd(wd_example)

#Import with import_files function that calls the list_csv function (filetype="csv")
viral_genome<-import_files(wd_example, filetype = "csv")

#Cleaning...
viral_vc<-cleaning(viral_genome, output_from="vcontact2")

#Creation of a abundance matrix with abundance_vc function
VC_mx_abs<-MetaViral::abundance_vc(viral_vc=viral_vc, abuntype = "absolute", taxa = "Genome",output_type = "matrix" )
VC_mx_rel<-MetaViral::abundance_vc(viral_vc=viral_vc, abuntype = "relative", taxa = "Genome",output_type = "matrix" )

#bio_statistics function
VC_stats<-bio_statistics(VC_mx_abs)
VC_stats<-bio_statistics(VC_mx_rel)
#-
#DRAMv
viral_annotations<-import_files(wd_example, filetype = "tsv")
viral_annotations<-cleaning(viral_annotations, output_from="DRAMv")
kegg<-MetaViral::db_exploring(viral_annotations = viral_annotations, database = "kegg")
peptidase<-MetaViral::db_exploring(viral_annotations = viral_annotations, database = "peptidase")
vog<-MetaViral::db_exploring(viral_annotations = viral_annotations, database = "vog")
all_db<-MetaViral::db_exploring(viral_annotations = viral_annotations, database = "all")

