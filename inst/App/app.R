library(shiny)
library(DT)
library(shinyFiles)
library(shinyjs)
library(MetaViral)
library(tibble)

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
        h4("Choose the directory with your data"),

        shinyFiles::shinyDirButton(id = "dataset", label = "Select path", title="Select path"),
        #shiny::selectInput(inputId = "dataset_type", label = "File type", choices = c("DRAM-v", "vConTACT2")),
        br(),
        br(),
        shiny::radioButtons(inputId = "dataset_type", label = "File type", choices = c("DRAM-v (tsv)", "vConTACT2 (csv)")),
        shiny::actionButton(inputId = "dataset_import", label = "Upload", icon = icon("log-out", lib = "glyphicon")),
        br(),
        hr(style = "border-color: black"),
        h5("...or use the example data"),
        shiny::actionButton(inputId = "ex_dram", label = "Example data (DRAM-v)"),
        shiny::actionButton(inputId = "ex_vcontact", label = "Example data (vConTACT2)"),
        br(),
        hr(style = "border-color: black"),
        shinyjs::hidden(shiny::selectInput(inputId = "amg_flags", label = "Flags to remove", choices = c("V", "M", "A", "P", "E", "K", "T", "F", "B"), multiple = TRUE)),
        shiny::actionButton(inputId = "dataset_clean", label = "Clean", icon = icon("clean", lib = "glyphicon")),
        br(),
        hr(style = "border-color: black"),
        br(),
        shinyjs::hidden(shiny::radioButtons(inputId = "dataset_matrix", label = "Abundance type to generate", choices = c("absolute", "relative"))),
        shinyjs::hidden(shiny::selectInput(inputId = "taxa_selection", label = "Taxon", choices = c("Order", "Family", "Genus", "Genome"))),
        shinyjs::hidden(shiny::actionButton(inputId = "dataset_abundance", label = "Abundance matrix", icon = icon("thumbnails-small", lib = "glyphicon"))),
        br(),
        br(),
        shinyjs::hidden(shiny::actionButton(inputId = "bio_stats", label = "Determine Biodiversity Indices",
                                            icon = icon("stats-circle", lib="glyphicon"))),

        shinyjs::hidden(shiny::selectInput(inputId = "database", label = "Database exploration",
                              choices=c("kegg", "pfam", "vog", "refseq", "cazy", "peptidase", "all"))),

        shinyjs::hidden(shiny::actionButton(inputId = "database_explore", label = "Explore"))



      ), #end of sidebarPanel
      shiny::mainPanel(fluidRow(
        h2("Dataset"),
        #img(src='logomv.png', align = "center"),

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


server <- function(input, output, session) {
  #Initiate variables
  typefile<-reactiveValues(df=NULL)
  data<-reactiveValues(df=NULL)
  shinyDirChoose(input, 'dataset', filetypes=c("csv", "tsv"), roots=c(wd = normalizePath("~")))
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
    cat(c("No", "directory", "selected"), sep = " ")
  })

  #Initiate observeEvents
  observeEvent(input$dataset_import, { #button to upload data
    req(input$dataset)

    if (input$dataset_type=="vConTACT2 (csv)"){
      typefile$df<-"csv"
    } else { #DRAM-v filetype
      typefile$df<-"tsv"
      shinyjs::show("amg_flags")
    }
    x<<-input$dataset
    wd<<-file.path(normalizePath("~"), paste(unlist(x$path[-1]), collapse = .Platform$file.sep))
    output$dir<-renderText({wd})

   data$df<-MetaViral::import_files(path=wd, filetype = typefile$df)
   output$imported_table<-DT::renderDataTable({
      data$df
    },
    options=(list(pageLength=10,scrollX=TRUE))
    )
   showNotification("Done")
   shinyjs::show("downloadRaw")
   shinyjs::show("downTypeRaw")

   printer$df<-data$df

  })#end observeEvent import

  observeEvent(input$ex_dram, { #button to upload example data (DRAM-v)
    typefile$df<-"tsv"
    data$df<-MetaViral::import_files(path="./www/ExampleData", filetype = "tsv")
    output$imported_table<-DT::renderDataTable({
      data$df
    },
    options=(list(pageLength=10,scrollX=TRUE))
    )
    showNotification("Done")
    shinyjs::show("downloadRaw")
    shinyjs::show("downTypeRaw")

    printer$df<-data$df
    shinyjs::show("amg_flags")

  })#end observeEvent example data (DRAM-v)

  observeEvent(input$ex_vcontact, { #button to upload example data (vcontact)
    typefile$df<-"csv"
    data$df<-MetaViral::import_files(path="./www/ExampleData", filetype = "csv")
    output$imported_table<-DT::renderDataTable({
      data$df
    },
    options=(list(pageLength=10,scrollX=TRUE))
    )
    showNotification("Done")
    shinyjs::show("downloadRaw")
    shinyjs::show("downTypeRaw")

    printer$df<-data$df

  })#end observeEvent example data (vcontact)

  observeEvent(input$dataset_clean, { #button to clean data
    tryCatch( #displays error on the app when crashing
      {
    if ( typefile$df == "csv"){
      output_from<-"vcontact2"
    }else{
      output_from<-"DRAMv"

    }
    data_clean$df<- MetaViral::cleaning(viral_data = data$df, output_from = output_from, remove_flags = input$amg_flags)
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
    abundance$df<-MetaViral::abundance_vc(viral_vc = data_clean$df,
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
    biostats<-MetaViral::bio_statistics(abundance$df)

    output$imported_table<-DT::renderDataTable({
      biostats
    },
    options=(list(pageLength=10,scrollX=TRUE))
    )
    showNotification("Done")
  })#end of observeEvent bio_stats

  observeEvent(input$database_explore, {
    data_base<-MetaViral::db_exploring(viral_annotations = data_clean$df, database = input$database )

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
