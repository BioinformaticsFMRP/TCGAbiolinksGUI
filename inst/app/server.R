# @title  Server side
# @description Server side - Download data from roadmap project
# @param input - input signal
# @param output - output signal
# @keywords internal
.biOMICsServer <- function(input, output,session) {

  .result <- reactiveValues(g_nbFiles = 0, g_downloaded = 0, df = NULL, downloading = 1)
  setwd(Sys.getenv("HOME"))
  rmapFolder   <- paste0(Sys.getenv("HOME"),"/GEO")
  encodeFolder <- paste0(Sys.getenv("HOME"),"/ENCODE")

  # ROADMAP TAB
  getRmapProject <- reactive ({ c(input$rmapProject,"Project","AND") })
  getRmapSearch  <- reactive ({ c(input$rmapSearch,NULL,NULL)        })
  getRmapType    <- reactive ({ c(input$rmapType,"Filter","AND")     })

  output$savedFiles <- renderUI({
    if(input$rmapDownloadBt){
      valueBoxOutput("rmapProgressBox", width = NULL)
    }
  })


    output$statusBox <- renderUI({
      invalidateLater(1000, session)
      valueBoxOutput("statusBoxText", width = NULL)
    })

  output$statusBoxText <- renderValueBox({
    if(.result$downloading == 2){
    valueBox(
      #h4(paste0("Downloaded  ",getNbFiles(), " files")), "Status", icon = icon("fa fa-spinner fa-spin"),
      h4(paste0("Downloaded files")), "Status", icon = icon("check-circle"),
      color = "green"
    )} else if (.result$downloading == 1) {
      valueBox(
        h4(paste0("Idle")), "Status", icon = icon("stop"),
        color = "red"
      )
    } else {
      valueBox(
        h4(paste0("Downloaded  ",getNbFiles(), " files")), "Status", icon = icon("fa fa-spinner fa-spin"),
        color = "yellow"
      )

    }
  })

  output$tbl <- DT::renderDataTable({
    if(input$rmapDownloadBt){  # trigger this function by pressing download button
      .result$downloading <- 3
      query <- eGeo(
                    isolate (getRmapSearch ()),
                    isolate (getRmapProject()),
                    isolate (getRmapType   ())
      )
      #print(query)
      geoDownloader(query,rmapFolder)
      if(!is.null(.result$df)) DT::datatable(.result$df)
    }

  })

  # ENCODE TAB
  getEncodeAssay    <- reactive({  input$assay    })
  getEncodeTarget   <- reactive({  input$target   })
  getEncodeFType    <- reactive({  input$ftype    })
  getEncodeSample   <- reactive({  input$sample   })
  getEncodeAssembly <- reactive({  input$assembly })
  getEncodeSearch   <- reactive({  input$search   })
  getNbFiles        <- reactive({  .result$g_nbFiles  })


  output$tblEncode <- DT::renderDataTable({
    if(input$encodeDownloadBt){  # trigger this function by pressing download button
      encodeDownloader(isolate(getEncodeSearch()),
                       "experiment",
                       isolate(getEncodeTarget()),
                       isolate(getEncodeSample()),
                       isolate(getEncodeFType()),
                       isolate(getEncodeAssay()),
                       isolate(getEncodeAssembly()),
                       encodeFolder
      )
      if(!is.null(.result$df)) DT::datatable(.result$df)
    }
  })


  output$savedPath <- renderValueBox({
    valueBox(
      h4(paste0(encodeFolder), "Output directory", icon = icon("folder-open"),
      color = "blue", width = 4
    )
  })

  output$savedPath2 <- renderValueBox({
    valueBox(
      h4(paste0(rmapFolder), "Output directory", icon = icon("folder-open"),
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



  # TCGA
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
          need(!is.null(tcga), "Please upload file with TCGA bar codes")
          )
        validate(
          need(!is.null(filter), "Please upload file with filtered bar codes")
        )

        barCode  <- .getBarCode (tcga, filter)
        fileConn <- file(filename)
        writeLines (barCode, fileConn)
        close (fileConn)
      }
    )

}
