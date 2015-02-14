#' @title  Server side
#' @description Server side - Download data from roadmap project
#' @name geoServer
#' @keywords internal
#' @param input - input signal
#' @param output - output signal
biOMICsServer <- function(input, output) {

  getProject <- reactive ({ c(input$project,"Project","AND") })
  getSearch  <- reactive ({ c(input$search,NULL,NULL)        })
  getType    <- reactive ({ c(input$type,"Filter","AND")     })

  output$value <- renderPrint({
    if(input$downloadBt){  # trigger this function by pressing download button
      query <- eGeo(
                    isolate (getSearch ()),
                    isolate (getProject()),
                    isolate (getType   ())
      )
      print(query)
      geoDownloader(query,"../download")
      print(paste0("Found: ",nbFiles," files."))
      print(paste0("Downloaded: ",nbFilesDownloaded," files."))
      #geoDownloader(query,"../download")
    }else{
      print("")
    }

  })
}
