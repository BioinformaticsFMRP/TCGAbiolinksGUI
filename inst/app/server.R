library(shiny)
library(shinyFiles)
#' @title  Server side
#' @description Server side - Download data from roadmap project
#' @param input - input signal
#' @param output - output signal
#' @keywords internal
biOMICsServer <- function(input, output, session) {
    setwd(Sys.getenv("HOME"))
    volumes <- c('Working directory'=getwd())
    volumes <- c('R Installation'=R.home())
    shinyDirChoose(input, 'directory', roots=volumes, session=session, restrictions=system.file(package='base'))
    output$directorypath <- renderPrint({parseDirPath(volumes, input$directory)})

    #----------------- Ontology
    getontProject <- reactive({c(input$ontProject,"Project","AND") })
    getontSearch  <- reactive({c(input$ontSearch,NULL,NULL)        })
    getontType    <- reactive({c(input$ontType,"Filter","AND")     })

    output$savedFiles <- renderUI({
        if (input$ontSearchDownloadBt ) {
            valueBoxOutput("ontProgressBox", width = NULL)
        }
    })

    output$ontSearchLink <-  renderText({
        if (input$ontSearchDownloadBt ) {
            link <- c()
            print(input$allRows)
            df <- data.frame(database = input$allRows[seq(1, length(input$allRows), 4)],
                             ID = input$allRows[seq(2, length(input$allRows), 4)],
                             Sample = input$allRows[seq(3, length(input$allRows), 4)],
                             Experiment = input$allRows[seq(4, length(input$allRows), 4)])
            print(df)
            biOmicsDownload(df)

        }}
    )

    output$ontSearchtbl <- renderDataTable({
        indexes <- c()
        query <- data.frame()
        if (input$ontSearchBt) {
            link <- c()
            term <- input$ontSamplesFilter
            exp <-  input$ontExpFilter
            print(term)
            print(exp)
            query <- biOmicsSearch(input$ontSamplesFilter,
                                   experiment = input$ontExpFilter)
            # improve using subset - subset(data,selection,projection)
        }
        query
    },
    options = list(pageLength = 10,
                   scrollX = TRUE,
                   jQueryUI = TRUE,
                   pagingType = "full",
                   lengthMenu = c(10, 50, 100, -1),
                   language.emptyTable = "No results found",
                   "dom" = 'T<"clear">lfrtip',
                   "oTableTools" = list(
                       "sSelectedClass" = "selected",
                       "sRowSelect" = "os",
                       "sSwfPath" = paste0("//cdnjs.cloudflare.com/ajax/",
                                           "libs/datatables-tabletools/",
                                           "2.2.3/swf/copy_csv_xls.swf"),
                       "aButtons" = list(
                           list("sExtends" = "collection",
                                "sButtonText" = "Save",
                                "aButtons" = c("csv","xls")
                           )
                       )
                   )
    ), callback = "function(table) {
      table.on('click.dt', 'tr', function() {
        Shiny.onInputChange('allRows',
                            table.rows('.selected').data().toArray());
      });

    }"
    )

}
