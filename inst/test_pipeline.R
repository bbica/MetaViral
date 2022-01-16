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


#DRAMv
viral_annotations<-import_files(wd_example, filetype = "tsv")
viral_annotations_c<-cleaning(viral_annotations, output_from="DRAMv", remove_flags=c("F", "T", "P", "A"), max_aux_score = 3)
#the default parameter for remove_flags is remove_flags=c("F", "T", "P", "A");
#use remove_flags=c() to perform a cleaning without filtering out any flag on the amg_flags collumn, in the dramv dataset(viral_annotations);
#input the maximum auxiliary score in max_aux_score (lower values correspond to a greater confidence of virus prediction); default set to 3, this will discard aux.scores of 4 and 5

#db_exploring
kegg<-MetaViral::db_exploring(viral_annotations = viral_annotations, database = "kegg") #this action uses the KEGGREST package, it will take a bit longer (https://www.kegg.jp/)
peptidase<-MetaViral::db_exploring(viral_annotations = viral_annotations, database = "peptidase") #peptidade refers to the MEROPS (https://www.ebi.ac.uk/merops/)
vog<-MetaViral::db_exploring(viral_annotations = viral_annotations, database = "vog") #http://vogdb.org/
pfam<-MetaViral::db_exploring(viral_annotations = viral_annotations, database = "pfam") #https://pfam.xfam.org/
refseq<-MetaViral::db_exploring(viral_annotations = viral_annotations, database = "refseq") #RefSeq viral (https://www.ncbi.nlm.nih.gov/genome/viruses/)
cazy<-MetaViral::db_exploring(viral_annotations = viral_annotations, database = "cazy") #dbCAN, automated carbohydrate-active enzyme (CAZyme) annotation (http://bcb.unl.edu/dbCAN2/)
all_db<-MetaViral::db_exploring(viral_annotations = viral_annotations, database = "all")
