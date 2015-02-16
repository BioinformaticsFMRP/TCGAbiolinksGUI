#' @title  Server side
#' @description Server side - Download data from roadmap project
#' @name geoServer
#' @keywords internal
#' @param input - input signal
#' @param output - output signal
biOMICsServer <- function(input, output) {

  source("globals.R")
  # ROADMAP TAB
  getRmapProject <- reactive ({ c(input$rmapProject,"Project","AND") })
  getRmapSearch  <- reactive ({ c(input$rmapSearch,NULL,NULL)        })
  getRmapType    <- reactive ({ c(input$rmapType,"Filter","AND")     })


  output$savedFiles <- renderUI({
    if(input$rmapDownloadBt)  valueBoxOutput("rmapProgressBox", width = NULL)
  })


  output$rmapReturn <- renderPrint({

    if(input$rmapDownloadBt){  # trigger this function by pressing download button
      query <- eGeo(
                    isolate (getRmapSearch ()),
                    isolate (getRmapProject()),
                    isolate (getRmapType   ())
      )
      #print(query)
      geoDownloader(query,"../download")
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
