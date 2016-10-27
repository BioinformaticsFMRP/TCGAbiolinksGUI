#' @title TCGAbiolinks GUI
#' @description Calls UI interface
#' @examples
#' \dontrun{
#'    TCGAbiolinksGUI()
#' }
#' @name TCGAbiolinksGUI
#' @import shiny shinyFiles shinydashboard downloader
#' TCGAbiolinks ggplot2 shinyBS stringr ggrepel pathview ELMER grid
#' clusterProfiler parallel readr data.table
#' @importFrom SummarizedExperiment SummarizedExperiment values rowRanges colData<- assay colData
#' @importFrom  colourpicker colourInput
#' @importFrom shinyjs hide show toggle useShinyjs
#' @export
#' @return Open a connection to shiny
TCGAbiolinksGUI <- function() {
    shiny::runApp(system.file("app", package = "TCGAbiolinksGUI"),launch.browser=TRUE)
}
