library(shiny)

ui <- navbarPage(
    source(file.path("ui", "ui.R"),  local = TRUE)$value
)

server <- function(input, output, session) {
    # Include the logic (server) for each tab
    source(file.path("server", "server.R"),  local = TRUE)$value
}
shinyApp(ui = ui, server = server)
