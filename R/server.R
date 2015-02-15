#' @title  Server side
#' @description Server side - Download data from roadmap project
#' @name geoServer
#' @keywords internal
#' @param input - input signal
#' @param output - output signal
biOMICsServer <- function(input, output) {

  # ROADMAP TAB
  getRmapProject <- reactive ({ c(input$rmapProject,"Project","AND") })
  getRmapSearch  <- reactive ({ c(input$rmapSearch,NULL,NULL)        })
  getRmapType    <- reactive ({ c(input$rmapType,"Filter","AND")     })

  output$rmapReturn <- renderPrint({
    if(input$rmapDownloadBt){  # trigger this function by pressing download button
      query <- eGeo(
                    isolate (getRmapSearch ()),
                    isolate (getRmapProject()),
                    isolate (getRmapType   ())
      )
      print(query)
      geoDownloader(query,"../download")
      print(paste0("Found: ", nbFiles," files."))
      print(paste0("Downloaded: ", nbFilesDownloaded," files."))
    }else{
      print("")
    }
  })

  # ENCODE TAB
  getEncodeAssay    <- reactive({  input$assay    })
  getEncodeTarget   <- reactive({  input$target   })
  getEncodeFType    <- reactive({  input$ftype    })
  getEncodeSample   <- reactive({  input$sample   })
  getEncodeAssembly <- reactive({  input$assembly })
  getEncodeSearch   <- reactive({  input$search   })


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
    }else{
      print("")
    }
  })

  output$rmapProgressBox <- renderValueBox({
    valueBox(
      paste0(0 + (100*nbFilesDownloaded)/nbFiles, "%"), "Progress", icon = icon("list"),
      color = "yellow"
    )
  })
}
