library(shiny)
library(shinyjs)
library(shinydashboard)
library(DT)
library(RSQLite)
library(DBI)
library(plyr)
library(ggplot2)
library(dplyr)
library(patchwork)
library(stringr)

dbPath <- "C:\\Users\\hemalab\\SQLiteStudio\\CovidDb_v2.db"
con <- dbConnect(SQLite(), dbPath)

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
    tabPanel(title = "Mutations",
             selectInput("protein", "Select Protein", ''),
             selectInput("mutation", "Select Mutation", '', multiple=TRUE),
             dateRangeInput("dates",label = "Select Date Range", start = '2020-01-01'),
             radioButtons("interval", "Select Time interval", choices = c("Month" = "m", "Week" = "W")),
             plotOutput("timelapse"),
             actionButton("test", "test stuff"),
             actionButton("submit", "Submit", style="float:right")
    ),
    tabPanel(title = "Download",
             h2("Download Fasta sequences and their respective annotation data"),
             fluidRow(
               column(width = 3, disabled(selectInput("protein2", "Select Protein", ''))),
               column(width = 3, disabled(selectInput("mutation2", "Select Mutation", '',multiple =TRUE))),
               column(width = 3, radioButtons("filterByMutation", "Filter samples by 1 Specific Mutation ", choices = c("Yes" = "Y", "No" = "N"), selected = "N"))),
             dateRangeInput("downloadDates",label = "Select Date Range", start = '2020-01-01'),
             radioButtons("filterFailed", "Filter QC failed samples", choices = c("Yes" = "Y", "No" = "N"), selected = "Y"),
             actionButton("download", "Download", style="float:right",),
             fluidRow(
               htmlOutput("googleVis")
             )
    )
  ),

)

server <- function(input, output, session) {

  countGenderQuery <- 'select Gender, count(Gender) from Samples where Gender is not null group by Gender'
  countVaccinatedQuery <- 'select count(IsVaccinated) from Samples where IsVaccinated = 1 group by IsVaccinated'
  ageRangeQuery <- 'select Min(Age), Max(Age) from Samples;'
  pangoLineageQuery <- 'select PangoLineage, count(PangoLineage) from Samples where PangoLineage is not null group by PangoLineage order by count(PangoLineage) desc'

  # A reactive expression that returns a whole intervals-occurences(with zeros)dataframe
  occurencesDf <- reactive({
    dateFormat <- paste0('%Y-%', input$interval)
    intervalQuery  <-  sqlInterpolate(con,
                                      'select strftime(?intervalUnit, SamplingDate) Interval from Samples where SamplingDate between ?dateLow and ?dateHigh group by Interval',
                                      intervalUnit = dateFormat ,dateLow = input$dates[1], dateHigh = input$dates[2])
    intervals <- getDfFromQuery(dbPath, intervalQuery)
    intervals['Occurences'] = 0
    intervals
  })

  # A reactive expression that returns a dataframe with time intervals and selected mutation occurences
  intervaldf <- eventReactive(list(input$dates[1], input$dates[2], input$interval), {
    dateFormat <- paste0('%Y-%', input$interval)
    intervalQuery  <-  sqlInterpolate(con,
                                      'select strftime(?intervalUnit, SamplingDate) Interval from Samples where SamplingDate between ?dateLow and ?dateHigh group by Interval',
                                      intervalUnit = dateFormat ,dateLow = input$dates[1], dateHigh = input$dates[2])

    intervals <- getDfFromQuery(dbPath, intervalQuery)
    intervals['Occurences'] = 0

    mutation <- paste0(input$protein, ":",input$mutation[1])
    input <- input$interval
    mutationTimelapseQuery <- sqlInterpolate(con, 'select strftime(?intervalUnit, SamplingDate) interval, count(MutationSample.MutationsMutationId)  from Samples
                                             inner join MutationSample on MutationSample.SamplesSampleId = Samples.SampleId
                                             where MutationSample.MutationsMutationId = ?mut
                                             group by interval',intervalUnit = dateFormat, mut = mutation)
    df <- getDfFromQuery(dbPath, mutationTimelapseQuery)
    for(i in intersect(intervals$Interval, df$interval)){
      intervals[intervals$Interval == i, 'Occurences']  = df[df$interval == i, 'count(MutationSample.MutationsMutationId)']
    }
    intervals
  })

  # A reactive expression that returns the intervals-total mutations dataframe
  totalMutationsPerIntervalDf <- eventReactive(input$submit,{
    dateFormat <- paste0('%Y-%', input$interval)
    #query <- sqlInterpolate(con,
    #                        'select strftime(?intervalUnit, SamplingDate) interval, count(SampleId) as TotalSamples from Samples
    #                         where Samples.SamplingDate between ?dateLow and ?dateHigh
    #                         group by interval', intervalUnit = dateFormat ,dateLow = input$dates[1], dateHigh = input$dates[2])
    query <- sqlInterpolate(con,
                            'select strftime(?, SamplingDate) interval, count(SampleId) as TotalSamples from Samples
                             where Samples.SamplingDate between ? and ?
                             group by interval', dateFormat, input$dates[1], input$dates[2])
    df <- getDfFromQuery(dbPath, query)
    intervals <- occurencesDf()
    for(i in intersect(intervals$Interval, df$interval)){
      intervals[intervals$Interval == i, 'Occurences']  = df[df$interval == i, 'TotalSamples']
    }
    intervals
  })

  # Sandbox button
  observeEvent(input$test, {
    print(intervaldf())
  })

  # a reactive expression that return all Proteins
  proteinNames <- reactive({
    proteinNamesQuery <- 'select distinct Protein from Mutations'
    i <- getDfFromQuery(dbPath, proteinNamesQuery)
    i$Protein
  })
  # a reactive expression that returns all mutations in the selected Protein
  mutationNames <- reactive({
    mutationNamesQuery <- sqlInterpolate(con, "select distinct AminoacidMutation from Mutations where Protein = ?p", p = input$protein)
    m <- getDfFromQuery(dbPath, mutationNamesQuery)
    m$AminoacidMutation
  })
  mutationNames2 <- reactive({
    con <- dbConnect(SQLite(), dbPath)
    mutationNamesQuery <- sqlInterpolate(con, "select distinct AminoacidMutation from Mutations where Protein = ?p", p = input$protein2)

    m <- dbGetQuery(con, mutationNamesQuery)
    dbDisconnect(con)
    m$AminoacidMutation
  })

  # An observe event that updated the proteins select box
  observe({
    updateSelectizeInput(session, "protein", choices=proteinNames(), selected = NULL, server = TRUE)
  })
  # An observe event that updated the selected protein mutations select box
  observe({
    updateSelectizeInput(session, "mutation", choices=mutationNames(), server = TRUE)
  })
 observe({
   updateSelectizeInput(session, "protein2", choices=proteinNames(), selected = NULL, server = TRUE)
 })
 observe({
   updateSelectizeInput(session, "mutation2", choices=mutationNames2(), server = TRUE)
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
    lowDate <- reactive(input$downloadDates[1])
    highDate <- reactive(input$downloadDates[2])

    if (input$filterFailed == 'Y'){
      sqlfailed <- "and Samples.QcStatus = 'passed_qc'"
    }else{
      sqlfailed <- ""
    }
  print(sqlfailed)
    if (input$filterByMutation == 'Y'){
      mutation <- paste0(input$protein2, ":",input$mutation2[1])
      sql <- paste0("select Fasta,SampleId,PangoLineage,Region,Wave,ScorpioCall,QcStatus,SamplingDate,Gender,Age,IsVaccinated,Outcome from Samples
    inner join MutationSample on MutationSample.SamplesSampleId = Samples.SampleId
    where MutationSample.MutationsMutationId = ?mut and (Samples.SamplingDate between ?low and ?high) ",sqlfailed)
      getSamplesforDownloadQuery <- sqlInterpolate(con, sql, mut = mutation, low = lowDate(), high = highDate())
      print(getSamplesforDownloadQuery)
      df <- dbGetQuery(con, getSamplesforDownloadQuery)
    }else {
      sql <- paste0("select Fasta,SampleId,PangoLineage,Region,Wave,ScorpioCall,QcStatus,SamplingDate,Gender,Age,IsVaccinated,Outcome from Samples
    where (Samples.SamplingDate between ?low and ?high) ",sqlfailed)
      getSamplesforDownloadQuery <- sqlInterpolate(con, sql ,low = lowDate(), high = highDate())
      print(getSamplesforDownloadQuery)
      df <- dbGetQuery(con, getSamplesforDownloadQuery)
    }
    print(count(df))
    output$DownloadCount <- renderText({
      nrow(df)
      })
  })



  ############## Plot Strains by time interval#######################################################
  observeEvent(input$submit, {
    output$timelapse <- renderPlot({
      totalMutationsDf <- totalMutationsPerIntervalDf()
      df <- intervaldf()
      df[,'TotalMutations'] <- totalMutationsDf[,'Occurences']
      df[,'freq'] <- df[,'Occurences'] / totalMutationsDf[,'Occurences']
      for (row in 1:nrow(df)){
        val <- paste0(as.character(df[row,'Occurences']),'/',as.character(df[row,'TotalMutations' ]))
        df[row,'labels'] <- val
      }
      print(df)
      ggplot2::ggplot(data = df, aes(x = Interval, label=labels, group =1)) +
        geom_line( aes(y=freq)) +
        #geom_line( aes(y=Occurences)) +
        geom_text(aes(x = Interval, y=freq))+
        ylim(0,1)
    })

    # Renders the graph of mutation occurences
    output$googleVis <- googleVis::renderGvis({
      protein <- input$protein
      genomegff3<-"data\\NC_045512.2_annot.gff3"
      gff3<-read.delim(genomegff3,as.is=TRUE,skip=2,header=FALSE)
      plen <- gff3[gff3$V9 == protein, 5] - gff3[gff3$V9 == protein, 4]


      q <- "select AminoacidMutation as aamut, VariantClass as mutationtype, count(SampleId) as occurences
          from Samples
          inner join MutationSample on MutationSample.SamplesSampleId = Samples.SampleId
          inner join Mutations on  Mutations.MutationId  = MutationSample.MutationsMutationId
          where Protein='S' and SamplingDate between '2020-12-23' and '2021-06-08'
          group by AminoacidMutation order by occurences desc limit 100"
      df <- getDfForMutationFrequencyPLot('S',dbPath, q)
      df$occurences <- as.numeric(df$occurences)
      df$aa <- as.numeric(df$aa)
      maxValue <- max(df$occurences)
      print(df)
      Sys.sleep(0.3)
      gvisBubbleChart(df, idvar="aamut", xvar="aa", yvar="occurences", colorvar="mutationtype",
                      options = list(
                        title=paste0("Mutation frequency for protein ",protein," in user-provided dataset"),
                        bubble='{textStyle: {fontSize: 11, color: "black", bold: true}}',
                        colors='["red","cornflowerblue","green", "purple","#eba56c", "brown", "yellow", "#69ffcf"]',
                        sizeAxis='{maxSize: 5, maxValue: 100}',
                        height=600
                      ))
    })
  })


  # Renders the text for the total number of samples in the database
  output$genderCount <- renderText({
    countSamplesQuery <- "select count(*) from Samples"
    countSamples <- getDfFromQuery(dbPath, countSamplesQuery)
    countSamples$`count(*)`
  })
  # Renders the text for the number of samples that have passed the qc in the database
  output$passedCount <- renderText({
    countPassedQuery <- 'select count(*) from Samples where QcStatus = "passed_qc"'
    countPassed <- getDfFromQuery(dbPath, countPassedQuery)
    countPassed$`count(*)`
  })
  # Renders the text for the number of samples that have failed the qc in the database
  output$failedCount <- renderText({
    countFailedQuery <- 'select count(*) from Samples where QcStatus = "fail"'
    countFailed <- getDfFromQuery(dbPath, countFailedQuery)
    countFailed$`count(*)`
  })
  # Renders the table for the Regions
  output$regionTable = DT::renderDataTable({
    countRegionsQuery <- 'select Region, count(Region) from Samples where Region is not null group by Region order by count(Region) desc'
    countRegions <- getDfFromQuery(dbPath, countRegionsQuery)
    DT::datatable(countRegions)
  })
  # Renders the table for the scorpio call information
  output$strainsTable <- DT::renderDataTable({
    scoprioCallQuery <- 'select ScorpioCall, count(ScorpioCall) from Samples where ScorpioCall is not null group by ScorpioCall order by count(ScorpioCall) desc'
    strains <- getDfFromQuery(dbPath, scoprioCallQuery)
    DT::datatable(strains)
  })

}

getDfFromQuery <- function(databasePath, query) {
  con <- dbConnect(SQLite(), databasePath)
  df <- dbGetQuery(con, query)
  dbDisconnect(con)
  return (df)
}

getDfForMutationFrequencyPLot <- function(protein,databasePath, query){

  mutationFrequencyDf <- getDfFromQuery(databasePath, query)
  mutationFrequencyDf$aa = stringr::str_extract(mutationFrequencyDf$aamut, "\\d+")
  return(mutationFrequencyDf)
}

shinyApp(ui = ui, server = server)
