#' @title biofunctionfind
#' @description groups rows into a new column, a number of strings is used to find matching rows and marks them with a tag
#' #'
#' @param alldb dataframe
#' @import dplyr
#' @examples
#' \dontrun{
#' quiet(cat("test"))
#' quiet(warning("test"))
#' quiet(warning("test"), all=T)
#' }
#' # This is a function that suppresses log info
#' @export
biofunctionfind<-function(alldb){
  '%>%' <- tidyr::`%>%`
  #genome regulation
Genome<-c("dna", "rna","r_na", "ribonucle", "helicase", "cyclase", "helix", "gtp",
                       "translation", "recombina", "restriction", "transcript", "nucle", "regulator", "terminase",
                       "ribosomal", "clamp", "cytosine","strand", "expression", "dcmp")#Genome function, DNA related proteins

  DNA_related<-c("dna")
  RNA_related<-c("rna", "ribonucle", "r_na")
  Helicase<-c("helicase")
  Recombination<-c("recombina")
  Restriction<-c("restriction")
  Transcription<-c("transcript")
  Nucleal_activty<-c("nucle")
  Regulators<-c("regulator")
  GTP<-c("gtp")
  Translation<-c("translation")
  Ribossome<-c("ribosomal")
  Clamp<-("clamp")
  Cytosine<-c("cytosine")
  Strand_binding<-c("strand")
  Terminase<-c("terminase")
  Helix<-c("helix")
  Expression_factor<-c("expression")
  Regulator<-c("regulator")
  Polymerase<-c("polymerase")
  Primase<-c("primase")
  Endonuclease<-c("endonuclease", "endodeoxyribonuclease")
  Exonuclease<-c("exonuclease",  "exodeoxyribonuclease")
  dTTP_synthesis<-c("dcmp", "dutp", "udp", "dtmp", "dtdp" )


Hypothetical<-c("uncharacterized", "putative", "gene", "unknown", "unnamed", "hypothetical", "possible")


Energy<-c("cyclase", "atp", "nad", "gmp","phoh", "glyco", "gluco", "ose", "manno")

Cyclase<-c("cyclase")
ATP<-c("atp")
NAD<-c("nad")
GMP<-c("gmp")
phoH<-c("phoh")
Sugars<- c("glyco", "gluco", "ose", "manno")

Structure<-c("head", "neck", "plate", "capsid", "tail", "structu", "chaperon", "portal", "scaffold")

  Capsid<-c("head", "capsid")
  Neck<-c("neck")
  Baseplate<-c("plate")
  Tail<-c("tail")
  Other_structure<-c("structu")
  Chaperone<-c("chaperon")
  Portal<-c("portal")
  Scaffold<-c("scaffold")


Infection<-c("virulence", "virion", "lysozyme", "acetyl", "curli", "integrase", "spanin", "endolysin", "fusion", "spike","tape_measure")

  Virulence<-c("virulence")
  Virion<-c("virion")
  Curli<-c("curli" )
  Integrase<-c("integrase")
  Spanin<-c("spanin")#
  Endolysin<-c("endolysin")#used by dsDNA phages to disrupt their hosts and be released in the environment
  Lysozyme<- c("lysozyme", "acetyl")
  Fusion<-c("fusion")
  Spike<-c("spike")
  Tape_measure<-c("tape_measure")

Metals<-c("iron", "ferr", "zinc","cobalt", "copper")

Zinc<-c("zinc")

Cobalt<-c("cobalt")

Iron<-c("iron", "ferr")

Copper<-c("copper")


Other_function<-c("glyco", "gluco", "ose", "manno", "phosph", "synthase", "reductase", "peptidoglycan", "transpo", "transferase",
                  "chitin", "heat", "shock", "peptidase", "toxin")

Phosphate<-c("phosph")

Synthase<-c("synthase" )

Reductase<-c("reductase")

Peptidoglycan<-c("peptidoglycan")

Transport<-c("transpo", "transferase")

Chitin<-c("chitin")

Toxin<-c("toxin")

Heat_stability<-c("heat", "shock")

Peptidase<- c("peptidase")



alldb$major_function<-NA
alldb$minor_function<-NA
alldb$clean_hit<-janitor::make_clean_names(alldb$hit)
alldb<-alldb%>%
  matchrows( colname="clean_hit", newcol="major_function" , string=Other_function, tag="Other function")%>%
  matchrows( colname="clean_hit", newcol="major_function" , string=Metals, tag="Metals")%>%
  matchrows( colname="clean_hit", newcol="major_function" , string=Genome, tag="Genome")%>%
  matchrows( colname="clean_hit", newcol="major_function" , string=Structure, tag="Structure")%>%
  matchrows( colname="clean_hit", newcol="major_function" , string=Infection, tag="Infection")%>%
  matchrows( colname="clean_hit", newcol="major_function" , string=Energy, tag="Energy")


alldb<-alldb%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Sugars, tag="Sugars")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Transport, tag="Transport")%>%

  matchrows( colname="clean_hit", newcol="minor_function" , string=Chitin, tag="Chitin")%>%

  matchrows( colname="clean_hit", newcol="minor_function" , string=Peptidase, tag="Peptidase")%>%

  matchrows( colname="clean_hit", newcol="minor_function" , string=Toxin, tag="Toxin")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Phosphate, tag="Phosphate")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Peptidoglycan, tag="Peptidoglycan")%>%

  matchrows( colname="clean_hit", newcol="minor_function" , string=Heat_stability, tag="Heat-stability")%>%

  matchrows( colname="clean_hit", newcol="minor_function" , string=Reductase, tag="Reductase")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Synthase, tag="Synthase")%>%

  #Metals
  matchrows(colname="clean_hit", newcol="minor_function" , string=Iron, tag="Iron")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Zinc, tag="Zinc")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Cobalt, tag="Cobalt")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Copper, tag="Copper")%>%

  #Structure
  matchrows( colname="clean_hit", newcol="minor_function" , string=Capsid, tag="Capsid")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Neck, tag="Neck")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Baseplate, tag="Baseplate")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Tail, tag="Tail")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Other_structure, tag="Structure-related")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Portal, tag="Portal")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Chaperone, tag="Chaperone")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Scaffold, tag="Scaffold")%>%

  #Genome Regulation
  matchrows( colname="clean_hit", newcol="minor_function" , string=DNA_related, tag="DNA-related")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Helicase, tag="Helicase")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Recombination, tag="Recombination")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=GTP, tag="GTP")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Translation, tag="Translation")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Ribossome, tag="Ribossome")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Clamp, tag="Clamp")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Cytosine, tag="Cytosine")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Strand_binding, tag="Strand-binding")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Terminase, tag="Terminase")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Nucleal_activty, tag="Nucleal activity")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Helix, tag="Helix")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=RNA_related, tag="RNA-related")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Expression_factor, tag="Expression factor")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Restriction, tag="Restriction")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Transcription, tag="Transcription")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Regulator, tag="Regulator")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Endonuclease, tag="Endonuclease")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Exonuclease, tag="Exonuclease")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=dTTP_synthesis, tag="dTTP synthesis")%>%

  #Infection
  matchrows( colname="clean_hit", newcol="minor_function" , string=Virulence, tag="Virulence")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Virion, tag="Virion")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Lysozyme, tag="Lysozyme")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Curli, tag="Curli")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Integrase, tag="Integrase")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Spanin, tag="Spanin")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Endolysin, tag="Endolysin")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Fusion, tag="Fusion")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Spike, tag="Spike")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=Tape_measure, tag="Tape measure")%>%
  #Energy
  matchrows( colname="clean_hit", newcol="minor_function" , string=Cyclase, tag="Cyclase")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=ATP, tag="ATP")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=NAD, tag="NAD")%>%
  matchrows( colname="clean_hit", newcol="minor_function" , string=GMP, tag="GMP")%>%
  matchrows( colname="clean_hit", newcol="major_function" , string=Hypothetical, tag="Hypothetical")

  return(alldb)
}#end of function
