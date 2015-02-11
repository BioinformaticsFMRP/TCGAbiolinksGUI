library("shiny")

shinyServer(function(input, output) {
  
  getAssay    <- reactive({  input$assay    })
  getTarget   <- reactive({  input$target   })
  getFType    <- reactive({  input$ftype    })
  getSample   <- reactive({  input$sample   })
  getAssembly <- reactive({  input$assembly })
  getSearch   <- reactive({  input$search   })
  
  output$value <- renderPrint({ 
    if(input$downloadBt){ # trigger this function by pressing download button
      print("Downloading")  # Does it work?
      encodeDownload(isolate(getSearch()),
                     "experiment",
                     isolate(getTarget()),
                     isolate(getSample()), 
                     isolate(getFType()),
                     isolate(getAssay()),
                     isolate(getAssembly()),
                     "encodeDownload"
      )
      print("Downloaded")  
    }else{
      print("Not downloaded")
    }    
  })
  
})