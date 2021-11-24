library(shiny)
library(shinyjs)
library(shinydashboard)
library(DT)
library(RSQLite)
library(DBI)
library(plyr)
library(ggplot2)

dbPath <- "C:\\Users\\HemaLab\\Desktop\\PROJECTS\\shiny_virome\\CovidDb_v2.db"
con <- dbConnect(SQLite(), dbPath)
 
mutationsDf <- dbGetQuery(con, "select * from Mutations limit 5")
dbDisconnect(con)
#mutationsDf <- dbReadTable(con, "Mutations")
#mutationSampleDf <- dbReadTable(con, "MutationSample")
#dbDisconnect(conn = con)


ui <- fluidPage(useShinyjs(),
  navlistPanel(              
    tabPanel(title = "Home",
             fluidRow(
               column(width = 8, "No of total Samples :"),
               column(width = 4, textOutput(outputId = "genderCount"))
             ),
             fluidRow(
               column(width = 3,"QC passed:"),
               column(width = 3, textOutput(outputId = "passedCount")),
               column(width = 3, "QC failed:"),
               column(width = 3, textOutput(outputId = "failedCount"))
             ),
             DT::dataTableOutput("regionTable"),
             DT::dataTableOutput("strainsTable")
             
    ),
    #tabPanel(title = "Samples",
    #         DT::dataTableOutput(outputId = "regionTable")
    #         
    #),
    tabPanel(title = "Mutations",
             #uiOutput("proteins")
             selectInput("protein", "Select Protein", ''),
             selectInput("mutation", "Select Mutation", '', multiple=TRUE),
             dateRangeInput("dates",label = "Select Date Range", start = '2020-01-01'),
             radioButtons("interval", "Select Time interval", choices = c("Month" = "m", "Week" = "W")),
             plotOutput("timelapse"),
             actionButton("submit", "Submit", style="float:right")
    ),
    tabPanel(title = "Download",
             h2("Download Fasta sequences and their respective annotation data"),
             fluidRow(
               column(width = 3, disabled(selectInput("protein2", "Select Protein", ''))),
               column(width = 3, disabled(selectInput("mutation2", "Select Mutation", ''))),
               column(width = 3, radioButtons("filterByMutation", "Filter samples by Specific Mutation", choices = c("Yes" = "Y", "No" = "N"), selected = "N"))),
             dateRangeInput("downloadDates",label = "Select Date Range", start = '2020-01-01'),
             radioButtons("filterFailed", "Filter QC failed samples", choices = c("Yes" = "Y", "No" = "N"), selected = "Y"),
             actionButton("download", "Download", style="float:right",),
             fluidRow(
               column(width = 8, "No of Samples to download :"),
               column(width = 4, textOutput(outputId = "DownloadCount"))
             )
    )
  ),
  
)

server <- function(input, output, session) {
  countSamplesQuery <- "select count(*) from Samples"
  countPassedQuery <- 'select count(*) from Samples where QcStatus = "passed_qc"'
  countFailedQuery <- 'select count(*) from Samples where QcStatus = "fail"'
  countGenderQuery <- 'select Gender, count(Gender) from Samples where Gender is not null group by Gender'
  countVaccinatedQuery <- 'select count(IsVaccinated) from Samples where IsVaccinated = 1 group by IsVaccinated'
  countRegionsQuery <- 'select Region, count(Region) from Samples where Region is not null group by Region order by count(Region) desc'
  ageRangeQuery <- 'select Min(Age), Max(Age) from Samples;'
  pangoLineageQuery <- 'select PangoLineage, count(PangoLineage) from Samples where PangoLineage is not null group by PangoLineage order by count(PangoLineage) desc'
  scoprioCallQuery <- 'select ScorpioCall, count(ScorpioCall) from Samples where ScorpioCall is not null group by ScorpioCall order by count(ScorpioCall) desc'
  
  
  proteinNames <- reactive({
    proteinNamesQuery <- 'select distinct Protein from Mutations'
    con <- dbConnect(SQLite(), dbPath)
    i <- dbGetQuery(con, proteinNamesQuery)
    dbDisconnect(con)
    i$Protein
  })
  mutationNames <- reactive({
    con <- dbConnect(SQLite(), dbPath)
    mutationNamesQuery <- sqlInterpolate(con, "select distinct AminoacidMutation from Mutations where Protein = ?p", p = input$protein)
    
    m <- dbGetQuery(con, mutationNamesQuery)
    dbDisconnect(con)
    m$AminoacidMutation
  })
  
  observe({
    updateSelectizeInput(session, "protein", choices=proteinNames(), selected = NULL, server = TRUE)  
  })
  observe({
    updateSelectizeInput(session, "mutation", choices=mutationNames(), server = TRUE)
  })
  observe({
    updateSelectizeInput(session, "protein2", choices=proteinNames(), selected = NULL, server = TRUE)  
  })
  observe({
    updateSelectizeInput(session, "mutation2", choices=mutationNames(), server = TRUE)
  })
  ################Download Panel##########################################################################
  observeEvent(input$filterByMutation,{
    filterByMutation <- input$filterByMutation
    if (filterByMutation == 'Y'){
      enable("protein2")
      enable("mutation2")
    } else {
      disable("protein2")
      disable("mutation2")
    }
  })
  
  observeEvent(input$download,{
    on.exit(dbDisconnect(con), add = TRUE)
    con <- dbConnect(SQLite(), dbPath)
    lowDate <- input$downloadDates[1]
    highDate <- input$downloadDates[2]
    
    if (input$filterFailed == 'Y'){
      sqlfailed <- 'and QcStatus = passed_qc'
    }else{
      sqlfailed <- ''
    }
    
    if (input$filterByMutation == 'Y'){
      mutation <- paste0(input$protein2, ":",input$mutation2)
      getSamplesforDownloadQuery <- sqlInterpolate(con, "select Fasta,SampleId,PangoLineage,Region,Wave,ScorpioCall,QcStatus,SamplingDate,Gender,Age,IsVaccinated,Outcome from Samples
    inner join MutationSample on MutationSample.SamplesSampleId = Samples.SampleId
    where MutationSample.MutationsMutationId = ?mut and ?qc", mut = mutation, qc = sqlfailed) 
      df <- dbGetQuery(con, getSamplesforDownloadQuery)
    }else {
      getSamplesforDownloadQuery <- sqlInterpolate(con, "select Fasta,SampleId,PangoLineage,Region,Wave,ScorpioCall,QcStatus,SamplingDate,Gender,Age,IsVaccinated,Outcome from Samples
    where MutationSample.MutationsMutationId = ?mut ?qc", qc = sqlfailed) 
      df <- dbGetQuery(con, getSamplesforDownloadQuery)
    }
    
  })
  
  
  
  ############## Plot Strains by time interval#######################################################
  observeEvent(input$submit, {
    output$timelapse <- renderPlot({
      on.exit(dbDisconnect(con), add = TRUE)
      con <- dbConnect(SQLite(), dbPath)
      mutation <- paste0(input$protein, ":",input$mutation[1])
      input <- input$interval
      mutationTimelapseQuery <- sqlInterpolate(con, "select strftime('%Y-%m', SamplingDate) month, count(MutationSample.MutationsMutationId)  from Samples
    inner join MutationSample on MutationSample.SamplesSampleId = Samples.SampleId
    where MutationSample.MutationsMutationId = ?mut
    group by month", mut = mutation) 
      df <- dbGetQuery(con, mutationTimelapseQuery)
      ggplot2::ggplot(data = df, aes(x = month, y = `count(MutationSample.MutationsMutationId)`, group = 1)) + 
        geom_path() + 
        geom_point()
    })
  })
  #output$proteins <- renderUI({
  #  selectInput("protein", "Select Protein", choices = getProteinNames())
  #})
 ################### Home Tab calculations ################################################################
  output$genderCount <- renderText({
    on.exit(dbDisconnect(con), add = TRUE)
    con <- dbConnect(SQLite(), dbPath)
    countSamples <- dbGetQuery(con, countSamplesQuery)
    countSamples$`count(*)`
  })
  output$passedCount <- renderText({
    on.exit(dbDisconnect(con), add = TRUE)
    con <- dbConnect(SQLite(), dbPath)
    countPassed <- dbGetQuery(con, countPassedQuery)
    countPassed$`count(*)`
  })
  output$failedCount <- renderText({
    on.exit(dbDisconnect(con), add = TRUE)
    con <- dbConnect(SQLite(), dbPath)
    countFailed <- dbGetQuery(con, countFailedQuery)
    countFailed$`count(*)`
  })
  output$regionTable = DT::renderDataTable({
    on.exit(dbDisconnect(con), add = TRUE)
    con <- dbConnect(SQLite(), dbPath)
    countRegions <- dbGetQuery(con, countRegionsQuery)
    DT::datatable(countRegions)
  })
  output$strainsTable <- DT::renderDataTable({
    on.exit(dbDisconnect(con), add = TRUE)
    con <- dbConnect(SQLite(), dbPath)
    strains <- dbGetQuery(con, scoprioCallQuery)
    DT::datatable(strains)
  })
  
  
}

getProteinNames <- function(){
  proteinNamesQuery <- 'select distinct Protein from Mutations'
  con <- dbConnect(SQLite(), dbPath)
  i <- dbGetQuery(con, proteinNamesQuery)
  dbDisconnect(con)
  print(i) 
  return(i$Protein)
}


shinyApp(ui = ui, server = server)