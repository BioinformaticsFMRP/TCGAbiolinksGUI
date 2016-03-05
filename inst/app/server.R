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
    shinyDirChoose(input, 'folder', roots=volumes, session=session, restrictions=system.file(package='base'))
    output$directorypath <- renderText({parseDirPath(volumes, input$folder)})

    #----------------- Ontology

    output$savedFiles <- renderUI({
        if (input$ontSearchDownloadBt ) {
            valueBoxOutput("ontProgressBox", width = NULL)
        }
    })

    output$ontSearchLink <-  renderText({
        getPath <- parseDirPath(volumes, input$folder)
        getFtype <- input$ontftypeFilter
        if (input$ontSearchDownloadBt ) {
            link <- c()

            df <- data.frame(database = input$allRows[seq(1, length(input$allRows), 5)],
                             ID = input$allRows[seq(2, length(input$allRows), 5)],
                             Sample = input$allRows[seq(3, length(input$allRows), 5)],
                             Experiment = input$allRows[seq(4, length(input$allRows), 5)],
                             organism = input$allRows[seq(5, length(input$allRows), 5)])
            if(length(getPath) > 0) {
                biOmicsDownload(df, path = getPath, enc.file.type = getFtype, rmap.file.type = getFtype)
            } else {
                biOmicsDownload(df, enc.file.type = getFtype, rmap.file.type = getFtype)
            }
            print("End of download")
        }}
    )

    dataInput <- reactive({
        indexes <- c()
        query <- data.frame()
            link <- c()
            term <- input$ontSamplesFilter
            exp <-  input$ontExpFilter
            query <- biOmicsSearch(input$ontSamplesFilter,
                                   experiment = input$ontExpFilter)
            # improve using subset - subset(data,selection,projection)
            return(list(system = names(sort(table(biosample.encode[biosample.encode$biosample %in% query$Sample,]$system),decreasing = T)[1]),
                        result = query))
    })

    output$system <-  renderText({
        if (input$ontSearchBt) {
            paste0("System: ", dataInput()$system)

        }
    })
    output$ontSearchtbl <- renderDataTable({
        if (input$ontSearchBt) {
            dataInput()$result
        }
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
