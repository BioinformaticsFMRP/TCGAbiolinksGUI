shinyServer(function(input, output) {

  # You can access the values of the widget (as a vector)
  # with input$checkGroup, e.g.
  output$assay <- renderPrint({ input$assay })
  output$target <- renderPrint({ input$target })
  
})