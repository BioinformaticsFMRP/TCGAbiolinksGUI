#' @title  Server side 
#' @description Server side - Download data from ENCODDE project
#' @name encodeDownloaderServer
#' @keywords internal
#' @param input - input signal
#' @param output - output signal
encodeServer <- function(input, output) {
  
  getAssay    <- reactive({  input$assay    })
  getTarget   <- reactive({  input$target   })
  getFType    <- reactive({  input$ftype    })
  getSample   <- reactive({  input$sample   })
  getAssembly <- reactive({  input$assembly })
  getSearch   <- reactive({  input$search   })
  
  output$value <- renderPrint({ 
    if(input$downloadBt){  # trigger this function by pressing download button
      encodeDownloader(isolate(getSearch()),
                       "experiment",
                       isolate(getTarget()),
                       isolate(getSample()), 
                       isolate(getFType()),
                       isolate(getAssay()),
                       isolate(getAssembly()),
                       "../download"
      )
    }else{
      print("")
    }
    
  })
}
