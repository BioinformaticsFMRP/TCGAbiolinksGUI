#' @title TCGAbiolinks GUI
#' @description Calls UI interface
#' @param run Used to control the examples.
#' @examples
#' \dontrun{
#'    TCGAbiolinksGUI()
#' }
#' TCGAbiolinksGUI(run = FALSE)
#' @name TCGAbiolinksGUI
#' @import shiny shinyFiles shinydashboard downloader
#' TCGAbiolinks ggplot2 shinyBS stringr ggrepel pathview ELMER grid
#' clusterProfiler parallel readr data.table maftools
#' @importFrom SummarizedExperiment SummarizedExperiment values rowRanges colData<- assay colData
#' @importFrom colourpicker colourInput
#' @importFrom shinyjs hide show toggle useShinyjs
#' @importFrom plotly renderPlotly plot_ly layout config plotlyOutput
#' @export
#' @return Open a connection to shiny
TCGAbiolinksGUI <- function(run = TRUE) {
    if(run) suppressMessages(shiny::runApp(system.file("app", package = "TCGAbiolinksGUI"),launch.browser=TRUE))
}
