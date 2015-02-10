library(shiny)

shinyUI(
  fluidPage(
    titlePanel("Teste"),
    
    sidebarLayout(
      position = "left",
      sidebarPanel( textInput(inputId = "text", label = NULL),
                    actionButton("go", label = "Go!")),
      mainPanel(
        textOutput("text1")
      )
    )
  )
)