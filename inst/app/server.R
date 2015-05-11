#' @title  Server side
#' @description Server side - Download data from roadmap project
#' @param input - input signal
#' @param output - output signal
#' @keywords internal
biOMICsServer <- function(input, output, session) {

    result <- reactiveValues(g_nbFiles = 0, g_downloaded = 0, df = NULL,
                             downloading = 1)
    setwd(Sys.getenv("HOME"))
    rmapDir   <- paste0(Sys.getenv("HOME"),"/GEO")
    encodeFolder <- paste0(Sys.getenv("HOME"),"/ENCODE")

    # ------------------------- ROADMAP TAB -------------------------
    getRmapProject <- reactive ({ c(input$rmapProject,"Project","AND") })
    getRmapSearch  <- reactive ({ c(input$rmapSearch,NULL,NULL)        })
    getRmapType    <- reactive ({ c(input$rmapType,"Filter","AND")     })

    output$savedFiles <- renderUI({
        if(input$rmapSearchDownloadBt){
            valueBoxOutput("rmapProgressBox", width = NULL)
        }
    })
    observe({
        if(input$selectDir){
            rmapDir <<- tcltk::tk_choose.dir(getwd(),
                                             paste0("Choose a suitable folder",
                                                    " to save the files"))
        }
        output$rmapDir <- renderText({rmapDir})
    })

    output$savedPath <- renderValueBox({
        input$selectDir
        valueBox(
            h4(rmapDir), "Output directory", icon = icon("folder-open"),
            color = "blue"
        )
    })

    output$statusBox <- renderUI({
        invalidateLater(1000, session)
        valueBoxOutput("statusBoxText", width = NULL)
    })


    output$statusBoxText <- renderValueBox({
        if(result$downloading == 2){
            valueBox(
                h4(paste0("Downloaded files")), "Status",
                icon = icon("check-circle"),
                color = "green"
            )} else if (result$downloading == 1) {
                valueBox(
                    h4(paste0("Idle")), "Status", icon = icon("stop"),
                    color = "red"
                )
            } else {
                valueBox(
                    h4(paste0("Downloaded  ",getNbFiles(), " files")),
                    "Status", icon = icon("fa fa-spinner fa-spin"),
                    color = "yellow"
                )

            }
    })

    output$rmapSearchLink <-  renderText({
        if(input$rmapSearchDownloadBt){
            link <- c()
            accession <- unlist(input$allRows)
            for (i in seq(1,length(accession), by = 6)){
                index <- which(roadmap.db$X..GEO.Accession ==  accession[i])
                link <- c(link,as.character(roadmap.db$GEO.FTP[index]))
            }
            geoDownloaderLinks(link,rmapDir, TRUE)
        }}
    )

    output$rmapSearchtbl <- renderDataTable({
        indexes <- c()
        if (input$rmapSearchBt){
            link <- c()
            # improve using subset - subset(data,selection,projection)
            for (i in seq_along(roadmap.db$Experiment)){
                if (is.null(input$rmapExpFilter)
                    || roadmap.db$Experiment[i]  %in% input$rmapExpFilter){
                    if (is.null(input$rmapSamplesFilter) ||
                            roadmap.db$Sample.Name[i] %in%
                            input$rmapSamplesFilter){
                        link <- c(link,as.character(roadmap.db$GEO.FTP[i]))
                        indexes <- c(indexes,i)
                    }
                }
            }
        }
        roadmap.db[indexes,c(1:5,8)]
    },
    options = list(pageLength = 10,
                   scrollX = TRUE,
                   jQueryUI= TRUE,
                   pagingType = "full",
                   lengthMenu = c(10, 50, 100, -1),
                   language.emptyTable = "No results found",
                   "dom" = 'T<"clear">lfrtip',
                   "oTableTools" = list(
                       "sSelectedClass"= "selected",
                       "sRowSelect"= "os",
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

    # ------------------------- ENCODE TAB -------------------------
    getEncodeAssay    <- reactive({  input$assay    })
    getEncodeTarget   <- reactive({  input$target   })
    getEncodeFType    <- reactive({  input$ftype    })
    getEncodeSample   <- reactive({  input$sample   })
    getEncodeAssembly <- reactive({  input$assembly })
    getEncodeSearch   <- reactive({  input$search   })
    getNbFiles        <- reactive({  result$g_nbFiles  })

    output$tblEncode <- renderDataTable({
        # trigger this function by pressing download button
        if(input$encodeDownloadBt){
            encodeDownloader(isolate(getEncodeSearch()),
                             "experiment",
                             isolate(getEncodeTarget()),
                             isolate(getEncodeSample()),
                             isolate(getEncodeFType()),
                             isolate(getEncodeAssay()),
                             isolate(getEncodeAssembly()),
                             encodeFolder
            )
            if(!is.null(result$df)) result$df
        }
    })



    output$savedPath2 <- renderValueBox({
        valueBox(
            h4(encodeFolder), "Output directory", icon = icon("folder-open"),
            color = "blue"
        )
    })


    output$rmapProgressBox <- renderValueBox({
        color <- "green"
        if(getNbFiles() == 0) color <- "red"

        valueBox(
            paste0(getNbFiles()), "Files saved", icon = icon("cloud-download"),
            color = color, width = 4
        )
    })

    # ------------------------- TCGA ----------------------
    tcgaBarCode  <- reactive({
        inFile <- input$file1
        if (is.null(inFile))
            return(NULL)
        read.csv(inFile$datapath, header = TRUE)
    })

    filterTcgaBarCode  <- reactive({
        inFile <- input$file2
        if (is.null(inFile))
            return(NULL)
        read.table(inFile$datapath, header = FALSE)
    })

    output$getTcgaBarCode <- downloadHandler(

        # This function returns a string which tells the client
        # browser what name to use when saving the file.
        filename = "barCode.txt",
        # This function should write data to a file given to it by
        # the argument 'file'.
        content = function(filename) {
            tcga <- tcgaBarCode()
            filter <- filterTcgaBarCode()
            validate(
                need(!is.null(tcga),
                     "Please upload file with TCGA bar codes")
            )
            validate(
                need(!is.null(filter),
                     "Please upload file with filtered bar codes")
            )

            barCode  <- getBarCode (tcga, filter)
            fileConn <- file(filename)
            writeLines (barCode, fileConn)
            close (fileConn)
        }
    )

}
