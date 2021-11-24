library(shiny)
library(shinyjs)
library(shinydashboard)
library(DT)
library(RSQLite)
library(DBI)
library(plyr)
library(ggplot2)

dbPath <- "C:\\Users\\hemalab\\SQLiteStudio\\CovidDb_v2.db"
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
             actionButton("test", "test stuff"),
             actionButton("submit", "Submit", style="float:right")
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
  
  intervaldf <- eventReactive(list(input$dates[1], input$dates[2], input$interval), {
    con <- dbConnect(SQLite(), dbPath)
    on.exit(dbDisconnect(con))
    dateFormat <- paste0('%Y-%', input$interval)
    intervalQuery  <-  sqlInterpolate(con,
                                      'select strftime(?intervalUnit, SamplingDate) Interval from Samples where SamplingDate between ?dateLow and ?dateHigh group by Interval',
                                      intervalUnit = dateFormat ,dateLow = input$dates[1], dateHigh = input$dates[2])
    print(intervalQuery)
    
    intervals <- dbGetQuery(con, intervalQuery)
    if (is.null(intervals)){
      
    }
    intervals['Occurences'] = 0
    
    mutation <- paste0(input$protein, ":",input$mutation[1])
    input <- input$interval
    mutationTimelapseQuery <- sqlInterpolate(con, 'select strftime(?intervalUnit, SamplingDate) interval, count(MutationSample.MutationsMutationId)  from Samples
    inner join MutationSample on MutationSample.SamplesSampleId = Samples.SampleId
    where MutationSample.MutationsMutationId = ?mut
    group by interval',intervalUnit = dateFormat, mut = mutation) 
    df <- dbGetQuery(con, mutationTimelapseQuery)
    for(i in intersect(intervals$Interval, df$interval)){
      
      intervals[intervals$Interval == i, 'Occurences']  = df[df$interval == i, 'count(MutationSample.MutationsMutationId)'] 
    }
    intervals
  })
  
  # Sandbox button
  observeEvent(input$test, {

    print(intervaldf())
  
  })
  
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
  
  observeEvent(input$submit, {
    output$timelapse <- renderPlot({
      df <- intervaldf()
      ggplot2::ggplot(data = df, aes(x = Interval, y = `Occurences`, group = 1)) + 
        geom_path() + 
        geom_point()
    })
  })
  #output$proteins <- renderUI({
  #  selectInput("protein", "Select Protein", choices = getProteinNames())
  #})

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