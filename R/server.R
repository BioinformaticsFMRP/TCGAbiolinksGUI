#' @title  Server side
#' @description Server side - Download data from roadmap project
#' @name geoServer
#' @keywords internal
#' @param input - input signal
#' @param output - output signal
#' debugging options(shiny.error=browser)
biOMICsServer <- function(input, output,session) {

  source("globals.R")
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
    if(result$downloading == 2){
    valueBox(
      #h4(paste0("Downloaded  ",getNbFiles(), " files")), "Status", icon = icon("fa fa-spinner fa-spin"),
      h4(paste0("Downloaded files")), "Status", icon = icon("check-circle"),
      color = "green"
    )} else if (result$downloading == 1) {
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
      result$downloading <- 3
      query <- eGeo(
                    isolate (getRmapSearch ()),
                    isolate (getRmapProject()),
                    isolate (getRmapType   ())
      )
      #print(query)
      geoDownloader(query,"../download")
      if(!is.null(result$df)) DT::datatable(result$df)
    }

  })

  # ENCODE TAB
  getEncodeAssay    <- reactive({  input$assay    })
  getEncodeTarget   <- reactive({  input$target   })
  getEncodeFType    <- reactive({  input$ftype    })
  getEncodeSample   <- reactive({  input$sample   })
  getEncodeAssembly <- reactive({  input$assembly })
  getEncodeSearch   <- reactive({  input$search   })
  getNbFiles        <- reactive({  result$g_nbFiles  })


  output$value <- renderPrint({
    if(input$encodeDownloadBt){  # trigger this function by pressing download button
      encodeDownloader(isolate(getEncodeSearch()),
                       "experiment",
                       isolate(getEncodeTarget()),
                       isolate(getEncodeSample()),
                       isolate(getEncodeFType()),
                       isolate(getEncodeAssay()),
                       isolate(getEncodeAssembly()),
                       "../download"
      )
      getNbFiles()
    }
  })


  output$savedPath <- renderValueBox({
    valueBox(
      h4(paste0("../downloads/")), "Output directory", icon = icon("folder-open"),
      color = "blue", width = 4
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
}
