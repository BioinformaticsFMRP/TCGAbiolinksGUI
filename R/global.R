#' @title TCGAbiolinks GUI
#' @description Calls UI interface
#' @examples
#'    TCGAbiolinksGUI()
#' @name TCGAbiolinksGUI
#' @import shiny shinyFiles shinydashboard downloader
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

setOptionsProgressBar <- function(title, label) {
    # if(is.linux() || is.mac() ) opb <-
    # pbapply::pboptions(type='tk', title = title, label = label)
    # else if(is.windows()) opb <- pbapply::pboptions(type='win',
    # title = title, label = label) else
    opb <- pbapply::pboptions(type = "txt", char = "+", title = title,
        label = label)
}

#' Get object from inside the package
#' @keywords internal
#' @export
get.obj <- function(obj){
    return(get(obj))
}
