library(shiny)
library(DT)
library(shinyjs)
library(tibble)
library(dplyr)
library(KEGGREST)
library(readr)
library(purrr)
library(tidyr)
library(gtools)
library(shinyFiles)
library(data.table)

#Before deploying run: (this will enable the detection of the Bioconductor package)
#library(BiocManager)
#options(repos = BiocManager::repositories())

list_csv <- function(flnm) {
  '%>%' <- tidyr::`%>%`
  return_namefile <-
    readr::read_csv(flnm) %>%
    dplyr::mutate(filename = flnm)
  return(return_namefile)
}

list_tsv <- function(flnm) {
  '%>%' <- tidyr::`%>%`
  return_namefile <-
    readr::read_tsv(flnm) %>%
    dplyr::mutate(filename = flnm)
  return(return_namefile)
}

import_files <- function(path, filetype="tsv"){
  '%>%' <- tidyr::`%>%`

  if (filetype=="tsv"){
    export_files <-
      list.files(path=path, pattern = paste0("*.", filetype),
                 full.names = T) %>%
      purrr::map_df(~list_tsv(.))
  } else if (filetype=="csv"){
    export_files <-
      list.files(path=path, pattern = paste0("*.", filetype),
                 full.names = T) %>%
      purrr::map_df(~list_csv(.))
  } else {
    warning("Filetype must be 'csv' or 'tsv'")
  }

  return(export_files)
}

cleaning<-function(viral_data, output_from="vcontact2", remove_flags=c("F", "T", "P", "A"), max_aux_score=3, ...){
  '%>%' <- tidyr::`%>%`
  if (output_from=="vcontact2"){
    viral_genome <- tidyr::separate(viral_data, col=filename, into=c("Biome", "Category"), sep="-(?=[^-])", extra = "merge")
    viral_genome <- tidyr::separate(viral_genome, col=Biome, into=c("junk", "Biome"), sep="(?=[A-z0-9])", extra = "merge")
    viral_genome <- tidyr::separate(viral_genome, col=Category, into=c("junk", "Category"), sep="(?=[0-9])", remove = FALSE)
    viral_genome <- tidyr::separate(viral_genome, col=Category, into=c("Category", "junk"), sep="[-]", extra = "merge")
    viral_genome$junk <- NULL
    viral_genome$Category[is.na(viral_genome$Category)] <- 0 #This enables the latter manipulations even if Category info is absent

    viral_vc<-viral_genome %>%
      dplyr::group_by(`VC Subcluster`) %>%
      dplyr::arrange(Genome, Family, Order, Genus, Biome)

    viral_vc<-viral_vc[!is.na(viral_vc$VC), ]

    viral_vc <- tidyr::separate(viral_vc, col=VC, into=c("VC", "VC2"), sep="\\D", extra = "merge")

    viral_vc<-transform(viral_vc, VC = as.numeric(VC))
    viral_vc<-transform(viral_vc, VC2 = as.numeric(VC2))
    viral_vc <- viral_vc[order(viral_vc$VC, viral_vc$VC2, decreasing = FALSE), ]

    return(viral_vc)



  }else if (output_from=="DRAMv") {
    #extract only the name of the biome from the filename column
    DRAMv_with_sources <- tidyr::separate(viral_data, col=filename, into=c("Biome", "junk"), sep="_(?=[^_]+$)", extra = "merge")
    DRAMv_with_sources <- tidyr::separate(DRAMv_with_sources, col=Biome, into=c("junk2", "Biome"), sep="(?=[A-z0-9])", extra = "merge")
    DRAMv_with_sources$junk <- NULL
    DRAMv_with_sources$junk2 <- NULL

    DRAMv_with_sources2 <- DRAMv_with_sources%>%
      tidyr::separate(col=viral_hit, into=c("viral_function", "protein_and_vname"), sep="\\s", extra="merge")%>%
      tidyr::separate(col=protein_and_vname, into=c("viral_function", "viral_tax"), sep="(?=\\[)")%>% #separates the column by the fist  "["
      #kegg
      tidyr::separate(col=kegg_hit, into=c("kegg_hit", "kegg_hit2"), sep=";")%>%
      tidyr::separate(col=kegg_hit, into=c("kegg_hit", "junk"), sep="(?=[^ ]+$)", extra = "merge")%>%
      #vogdb
      tidyr::separate(col=vogdb_description, into=c("vogdb_description", "junk"), sep=";")%>%
      tidyr::separate(col=vogdb_description, into=c("junk2", "vogdb_description"), sep="(?=\\s)", extra = "merge")%>%
      #pfam
      tidyr::separate(col=pfam_hits, into=c("pfam_hits", "pfam_id"), sep="(?=\\[)")%>%
      #peptidase
      tidyr::separate(col=peptidase_hit, into=c("junk", "peptidase_hit"), sep="-", extra = "merge")%>%
      tidyr::separate(col=peptidase_hit, into=c("peptidase_hit", "junk2"), sep="(?<=\\))", extra = "merge")%>%
      tidyr::separate(col=peptidase_hit, into=c("peptidase_hit", "peptidase_tax"), sep="(?=\\()", extra = "drop")%>%
      tidyr::separate(col=peptidase_tax, into=c("junk", "peptidase_tax"), sep="(?=\\w)", extra = "merge")%>%
      #Cazy
      tidyr::separate(col=cazy_hits, into=c("cazy_hits", "junk"), sep="(?=\\()", extra = "drop")%>%
      tidyr::separate(col=cazy_hits, into=c("junk2", "cazy_hits"), sep="\\s", extra = "merge")

    DRAMv_with_sources2$junk <- NULL
    DRAMv_with_sources2$junk2 <- NULL

    #Extracts the taxonomic predition from the square brackets
    DRAMv_with_sources2$viral_tax = sub(".*\\[([^][]+)].*", "\\1", DRAMv_with_sources2$viral_tax)

    #removes rows without any predictions
    DRAMv_with_sources3 <- DRAMv_with_sources2 %>%
      dplyr::filter_at( dplyr::vars(kegg_hit,viral_function,pfam_hits,cazy_hits,vogdb_description,peptidase_hit),dplyr::any_vars(!is.na(.)))

    #removes rows with NA in the amg_flags
    #DRAMv_with_sources3 <- DRAMv_with_sources3[!is.na(DRAMv_with_sources3$amg_flags),]
    #filter for only auxiliary scores inferior to 4
    DRAMv_with_sources4<-dplyr::filter(DRAMv_with_sources3, auxiliary_score<(max_aux_score+1))
    #filter out flags related with viral function
    if (rlang::is_empty(remove_flags)==TRUE){

    }else{
      for (i in 1:max(nchar(DRAMv_with_sources4$amg_flags))) {
        prm<-gtools::permutations(n=length(remove_flags),r=i,v=remove_flags, repeats.allowed = TRUE)
        prm<-apply(prm, 1, function(x) paste0(x, collapse=''))
        DRAMv_with_sources4<-DRAMv_with_sources4[ ! DRAMv_with_sources4$amg_flags %in% prm, ]
      }
    }

    viral_annotations<-DRAMv_with_sources4
    return(viral_annotations)

  }else {
    warning("output_from must be vcontact2 or DRAMv")
  }
}#end of function


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

db_exploring<-function(viral_annotations, database="all"){
  '%>%' <- tidyr::`%>%`
  if (database=="peptidase"){
    peptidase<-dplyr::select(viral_annotations, peptidase_hit, peptidase_id, peptidase_family, peptidase_RBH, peptidase_eVal, Biome, fasta)
    peptidase<-na.omit(peptidase)
    peptidase<-peptidase %>%
      dplyr::add_count(peptidase_hit,name="peptidase_counts", Biome)
    return(peptidase)

  }else if (database=="vog"){
    vog<-dplyr::select(viral_annotations, vogdb_description, Biome, fasta)
    vog<-na.omit(vog)
    vog<-vog %>%
      dplyr::add_count(vogdb_description,name="vog_counts", Biome)
    return(vog)

  }else if (database=="kegg"){
    kegg<-dplyr::select(viral_annotations, Biome, kegg_hit, kegg_id, fasta)
    kegg<-na.omit(kegg)
    kegg<-kegg %>%
      dplyr::add_count(kegg_hit,name="kegg_counts", Biome)

    kegg$Pathway<-NA
    kegg$Pathway2<-NA
    kegg$Pathway3<-NA

    for (i in 1:nrow(kegg)) {
      protein_id<-c(kegg$kegg_id[i])
      query<-tryCatch(KEGGREST::keggGet(protein_id), error=function(e) NULL)

      if (length(query[[1]]$PATHWAY)==1) {
        kegg$Pathway[i]<-query[[1]]$PATHWAY

      } else if (length(query[[1]]$PATHWAY)==2) {
        kegg$Pathway[i]<-query[[1]]$PATHWAY[1]
        kegg$Pathway2[i]<-query[[1]]$PATHWAY[2]

      } else if (length(query[[1]]$PATHWAY)==3) {
        kegg$Pathway[i]<-query[[1]]$PATHWAY[1]
        kegg$Pathway2[i]<-query[[1]]$PATHWAY[2]
        kegg$Pathway3[i]<-query[[1]]$PATHWAY[3]
      }
    }#end of for loop

    return(kegg)

  }else if (database=="pfam"){
    pfam<-dplyr::select(viral_annotations, Biome, pfam_hits, pfam_id, fasta)
    pfam<-na.omit(pfam)
    pfam<-pfam %>%
      dplyr::add_count(pfam_hits,name="pfam_counts", Biome)
    return(pfam)

  }else if (database=="refseq"){
    refseq<-dplyr::select(viral_annotations, Biome, viral_function, viral_id, viral_tax, viral_identity,
                          viral_bitScore, viral_eVal, fasta)
    refseq<-na.omit(refseq)
    refseq<-refseq %>%
      dplyr::add_count(viral_function,name="viral_counts", Biome)
    return(refseq)

  }else if (database=="cazy"){ #dbCAN, automated carbohydrate-active enzyme (CAZyme) annotation (http://bcb.unl.edu/dbCAN2/)
    cazy<-dplyr::select(viral_annotations, Biome, cazy_hits, fasta)
    cazy<-na.omit(cazy)
    cazy<- cazy %>%
      dplyr::add_count(cazy_hits, name="cazy_counts", Biome)
    return(cazy)

  }else if (database=="all"){
    peptidase<-dplyr::select(viral_annotations, peptidase_hit, peptidase_id, peptidase_family, peptidase_RBH, peptidase_eVal, Biome)
    peptidase<-na.omit(peptidase)
    peptidase<-peptidase %>%
      dplyr::add_count(peptidase_hit,name="peptidase_counts", Biome)
    peptidase<-unique(peptidase[, c("peptidase_hit","Biome","peptidase_counts")])

    vog<-dplyr::select(viral_annotations, vogdb_description, Biome)
    vog<-na.omit(vog)
    vog<-vog %>%
      dplyr::add_count(vogdb_description,name="vog_counts", Biome)
    vog<-unique(vog[, c("vogdb_description","Biome","vog_counts")])

    kegg<-dplyr::select(viral_annotations, Biome, kegg_hit, kegg_id)
    kegg<-na.omit(kegg)
    kegg<-kegg %>%
      dplyr::add_count(kegg_hit,name="kegg_counts", Biome)
    kegg<-unique(kegg[, c("kegg_hit","kegg_id","Biome","kegg_counts")])

    kegg$Pathway<-NA
    kegg$Pathway2<-NA
    kegg$Pathway3<-NA

    for (i in 1:nrow(kegg)) {
      protein_id<-c(kegg$kegg_id[i])
      query<-tryCatch(KEGGREST::keggGet(protein_id), error=function(e) NULL)

      if (length(query[[1]]$PATHWAY)==1) {
        kegg$Pathway[i]<-query[[1]]$PATHWAY

      } else if (length(query[[1]]$PATHWAY)==2) {
        kegg$Pathway[i]<-query[[1]]$PATHWAY[1]
        kegg$Pathway2[i]<-query[[1]]$PATHWAY[2]

      } else if (length(query[[1]]$PATHWAY)==3) {
        kegg$Pathway[i]<-query[[1]]$PATHWAY[1]
        kegg$Pathway2[i]<-query[[1]]$PATHWAY[2]
        kegg$Pathway3[i]<-query[[1]]$PATHWAY[3]
      }
    }#end of for loop

    pfam<-dplyr::select(viral_annotations, Biome, pfam_hits, pfam_id)
    pfam<-na.omit(pfam)
    pfam<-pfam %>%
      dplyr::add_count(pfam_hits,name="pfam_counts", Biome)
    pfam<-unique(pfam[, c("pfam_hits","Biome","pfam_counts")])


    refseq<-dplyr::select(viral_annotations, Biome, viral_function, viral_id, viral_tax, viral_identity,
                          viral_bitScore, viral_eVal)
    refseq<-na.omit(refseq)
    refseq<-refseq %>%
      dplyr::add_count(viral_function,name="viral_counts", Biome)
    refseq<-unique(refseq[, c("viral_function","Biome","viral_counts")])

    cazy<-dplyr::select(viral_annotations, Biome, cazy_hits)
    cazy<-na.omit(cazy)
    cazy<- cazy %>%
      dplyr::add_count(cazy_hits, name="cazy_counts", Biome)
    cazy<-unique(cazy[, c("cazy_hits","Biome","cazy_counts")])

    #merge all db
    database_all<-Reduce(function(x,y) merge(x,y,by="Biome",all=TRUE) ,list(kegg, pfam, vog, refseq, peptidase, cazy))
    return(database_all)

  }else {
    warning("choose one of the following: kegg, pfam, vog, refseq, cazy, peptidase or all")
  }
}#end of function


# UI #
ui <- shinyUI(
  fluidPage(
    tags$head(

      # Note the wrapping of the string in HTML()
      tags$style(HTML("
      @import url('https://fonts.googleapis.com/css2?family=Rajdhani:wght@600&display=swap');
      body {
        background-color: lightgrey;
        color: black;
      }
      h2 {
        font-family: 'Rajdhani', sans-serif;
      }
      .shiny-input-container {
        color: #474747;
      }"))
    ),
    shinyjs::useShinyjs(),
    div(
      id="myapp", #Start of header
      shiny::titlePanel(title = div(
        shiny::splitLayout(
          cellWidths = NULL,
          h2("MetaViral app", align = "left",
             img(height = 55, width = 75,
                 src = "logomv.png",
                 class = "pull-center")),
          div(
            style = "position:absolute;top:1em; right:1em;",
            a(img(height = 50, width = 50, src="logo.png"), href="https://github.com/bbica/MetaViral") ,
          ))
      ))

    ),#end of div

    shiny::sidebarLayout(
      shiny::sidebarPanel(
        fileInput("dataset", h4("Choose file(s):"),
                  accept = c(".tsv",".csv"),multiple = TRUE),
        #shinyFiles::shinyDirButton(id = "dataset", label = "Select path", title="Select path"),
        #shiny::selectInput(inputId = "dataset_type", label = "File type", choices = c("DRAM-v", "vConTACT2")),
        shiny::radioButtons(inputId = "dataset_type", label = "File type", choices = c("DRAM-v (tsv)", "vConTACT2 (csv)")),
        shiny::actionButton(inputId = "dataset_import", label = "Upload", icon = icon("log-out", lib = "glyphicon")),
        br(),
        hr(style = "border-color: black"),
        h5("...or use the example data"),
        shiny::actionButton(inputId = "ex_dram", label = "Example data (DRAM-v)"),
        shiny::actionButton(inputId = "ex_vcontact", label = "Example data (vConTACT2)"),
        br(),
        hr(style = "border-color: black"),
        #shinyjs::hidden(shiny::selectInput(inputId = "aux_score", label = "Max aux score", choices = c("1", "2", "3", "4", "5"), multiple = FALSE)),
        shinyjs::hidden(shiny::numericInput(inputId = "aux_score", label = "Max aux score", 4, min = 1, max = 5),
        verbatimTextOutput("value")),
        shinyjs::hidden(shiny::selectInput(inputId = "amg_flags", label = "Flags to remove", choices = c("V", "M", "A", "P", "E", "K", "T", "F", "B"), multiple = TRUE)),
        shinyjs::hidden(shiny::actionButton(inputId = "dataset_clean", label = "Process data", icon = icon("clean", lib = "glyphicon"))),
        br(),
        br(),
        hr(style = "border-color: black"),
        br(),
        shinyjs::hidden(shiny::radioButtons(inputId = "dataset_matrix", label = "Abundance type to generate", choices = c("absolute", "relative"))),
        shinyjs::hidden(shiny::selectInput(inputId = "taxa_selection", label = "Taxon", choices = c("Order", "Family", "Genus", "Genome"))),
        shinyjs::hidden(shiny::actionButton(inputId = "dataset_abundance", label = "Abundance matrix", icon = icon("thumbnails-small", lib = "glyphicon"))),
        br(),
        shinyjs::hidden(shiny::actionButton(inputId = "bio_stats", label = "Determine Biodiversity Indices", icon = icon("stats-circle", lib="glyphicon"))),
        #shiny::selectizeInput(inputId = "database", label = "Database exploration",
        #choices=c("kegg", "pfam", "vog", "viraldb", "peptidase", "all"), multiple=T),
        shinyjs::hidden(shiny::selectInput(inputId = "database", label = "Database exploration",
                                           choices=c("kegg", "pfam", "vog", "refseq", "cazy", "peptidase", "all"))),
        shinyjs::hidden(shiny::actionButton(inputId = "database_explore", label = "Explore"))



      ), #end of sidebarPanel
      shiny::mainPanel(fluidRow(
        h2("Dataset"),

        div(
          style = "position:absolute;top:1em; right:1em;",
          shinyjs::hidden(selectInput(inputId ="downTypeRaw", label = "Filetype", choices=c("rds", "csv"))),
          shinyjs::hidden(downloadButton(outputId = "downloadRaw", label = "Download"))
        ),
        br(),
        br(),
        h4("Filepath"),
        verbatimTextOutput("dir", placeholder = TRUE),
        uiOutput(outputId = "output_dataset"),
        h4("Data"),
        DT::dataTableOutput(outputId = "imported_table")
      ))
    )#end of sidebarLayout
  )#end of fluidPage
)#end of shintUI

# server #
server <- function(input, output, session) {
  #Initiate variables
  typefile<-reactiveValues(df=NULL)
  data<-reactiveValues(df=NULL)
  data_clean<-reactiveValues(df=NULL)
  abundance<-reactiveValues(df=NULL)
  printer<-reactiveValues(df=NULL)
  #Initiate placeholders
  output$imported_table <- DT::renderDataTable(
    {
      showData <- tibble(FileNames = vector("character"))
      isolate(showData)
    },
    rownames = FALSE,
    options = (list(scrollX = TRUE))
  )
  output$dir <- renderPrint({
    cat(c("No", "file(s)", "selected"), sep = " ")
  })

  #Initiate observeEvents
  observeEvent(input$dataset_import, { #button to upload data
    req(input$dataset)

    if (input$dataset_type=="vConTACT2 (csv)"){
      typefile$df<-"csv"
      shinyjs::hide("aux_score")
      shinyjs::hide("amg_flags")
      shinyjs::hide("database")
      shinyjs::hide("database_explore")
    } else { #DRAM-v filetype
      typefile$df<-"tsv"
      shinyjs::show("aux_score")
      shinyjs::show("amg_flags")
      shinyjs::hide("dataset_matrix")
      shinyjs::hide("taxa_selection")
      shinyjs::hide("dataset_abundance")
      shinyjs::hide("bio_stats")
    }
    x<<-input$dataset
    wd<<-file.path(normalizePath("~"), paste(unlist(x$path[-1]), collapse = .Platform$file.sep))
    output$dir<-renderText({wd})

    data$df<- map2_df(input$dataset$name, input$dataset$datapath,
                    ~data.table::fread(.y)%>% mutate(filename = .x))
    output$imported_table<-DT::renderDataTable({
      data$df
    },
    options=(list(pageLength=10,scrollX=TRUE))
    )
    showNotification("Done")
    shinyjs::show("downloadRaw")
    shinyjs::show("downTypeRaw")
    shinyjs::show("dataset_clean")
    printer$df<-data$df

  })#end observeEvente import

  observeEvent(input$ex_dram, { #button to upload example data (DRAM-v)
    typefile$df<-"tsv"
    data$df<-import_files(path="./www", filetype = "tsv")
    data$df$filename<-gsub("./www/", "", data$df$filename)
    output$imported_table<-DT::renderDataTable({
      data$df
    },
    options=(list(pageLength=10,scrollX=TRUE))
    )
    showNotification("Done")
    shinyjs::show("downloadRaw")
    shinyjs::show("downTypeRaw")
    shinyjs::show("dataset_clean")
    printer$df<-data$df
    shinyjs::show("aux_score")
    shinyjs::show("amg_flags")
    shinyjs::hide("dataset_matrix")
    shinyjs::hide("taxa_selection")
    shinyjs::hide("dataset_abundance")
    shinyjs::hide("bio_stats")

  })#end observeEvent example data (DRAM-v)

  observeEvent(input$ex_vcontact, { #button to upload example data (vcontact)
    typefile$df<-"csv"
    data$df<-import_files(path="./www", filetype = "csv")
    data$df$filename<-gsub("./www/", "", data$df$filename)
    output$imported_table<-DT::renderDataTable({
      data$df
    },
    options=(list(pageLength=10,scrollX=TRUE))
    )
    showNotification("Done")
    shinyjs::hide("aux_score")
    shinyjs::hide("amg_flags")
    shinyjs::hide("database")
    shinyjs::hide("database_explore")
    shinyjs::show("downloadRaw")
    shinyjs::show("downTypeRaw")
    shinyjs::show("dataset_clean")
    printer$df<-data$df

  })#end observeEvent example data (vcontact)


  observeEvent(input$dataset_clean, { #button to clean data
    tryCatch( #displays error on the app when crashing
      {
        if ( typefile$df == "csv"){
          output_from<-"vcontact2"
        }else{
          output_from<-"DRAMv"
          output$value <- renderText({ input$aux_score })

        }

        data_clean$df<- cleaning(viral_data = data$df, output_from = output_from, remove_flags = input$amg_flags, max_aux_score = input$aux_score )
        output$imported_table<-DT::renderDataTable({
          data_clean$df
        },
        options=(list(pageLength=10,scrollX=TRUE))
        )
        showNotification("Done")
        if (output_from == "vcontact2"){
          shinyjs::show("dataset_matrix")
          shinyjs::show("taxa_selection")
          shinyjs::show("dataset_abundance")
        }else{
          shinyjs::show("database")
          shinyjs::show("database_explore")
        }
        printer$df<-data_clean$df

      },
      error = function(err) { #displays error on the app when crashing
        showNotification(paste0(err), type = "err")
        message(err)
      }
    )#end of tryCatch (error display)
  })#end observeEvent clean

  observeEvent(input$dataset_abundance, {
    showNotification("Wait patiently...")
    abundance$df<-abundance_vc(viral_vc = data_clean$df,
                                          taxa = input$taxa_selection,
                                          abuntype = input$dataset_matrix,
                                          output_type = "matrix")

    output$imported_table<-DT::renderDataTable({
      abundance$df
    },
    options=(list(pageLength=10,scrollX=TRUE))
    )
    showNotification("Done")

    if (input$dataset_matrix=="absolute"){
      shinyjs::show("bio_stats")
    }
    printer$df<- abundance$df

  })#end of observeEvent abundance
  observeEvent(input$bio_stats, {
    biostats<-bio_statistics(abundance$df)

    output$imported_table<-DT::renderDataTable({
      biostats
    },
    options=(list(pageLength=10,scrollX=TRUE))
    )
    showNotification("Done")
  })#end of observeEvent bio_stats

  observeEvent(input$database_explore, {
    data_base<-db_exploring(viral_annotations = data_clean$df, database = input$database )

    output$imported_table<-DT::renderDataTable({
      data_base
    },
    options=(list(pageLength=10,scrollX=TRUE))
    )
    showNotification("Done")
  })



  # Download Raw Data
  output$downloadRaw <- downloadHandler(
    filename = function() {
      today <- toupper(format(Sys.Date(), format = "%d%b%y"))
      paste("ViralData_", today, ".", input$downTypeRaw, sep = "")
    },
    content = function(file) {
      tmp <- printer$df

      if (input$downTypeRaw == "rds") {
        saveRDS(object = tmp, file = file)
      } else {

        write.csv(x = tmp, file = file, quote = F, row.names = F, col.names = T)
      }
    }
  )

}
shinyApp(ui, server)
