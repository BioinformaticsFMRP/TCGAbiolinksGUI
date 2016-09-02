#' @title TCGAbiolinks GUI
#' @description Calls UI interface
#' @examples
#'    TCGAbiolinksGUI()
#' @name TCGAbiolinksGUI
#' @import shiny shinyFiles shinydashboard downloader colourpicker
#' TCGAbiolinks UpSetR ggplot2 shinyBS stringr ggrepel pathview ELMER grid
#' clusterProfiler parallel
#' @importFrom SummarizedExperiment SummarizedExperiment values rowRanges colData<- assay colData
#' @importFrom shinyjs hide show toggle useShinyjs
#' @export
#' @return Open a connection to shiny
TCGAbiolinksGUI <- function() {
    shiny::runApp(system.file("app", package = "TCGAbiolinksGUI"),launch.browser=TRUE)
}


is.windows <- function() {
    Sys.info()["sysname"] == "Windows"
}

is.mac <- function() {
    Sys.info()["sysname"] == "Darwin"
}

is.linux <- function() {
    Sys.info()["sysname"] == "Linux"
}

#' Get object from inside the package
#' @keywords internal
#' @export
#' @return Return object from package
#' @examples
#' maf.files <- get.obj("maf.files")
get.obj <- function(obj){
    if(is.null(obj)) return(NULL)
    return(get(obj))
}
