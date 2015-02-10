library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Hello bOMICs!"),
  
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      textInput("text", "Texto: ", "digite aqui")
      
    )
    
  )
  mainPanel(
    textOutput("text")
  )
))