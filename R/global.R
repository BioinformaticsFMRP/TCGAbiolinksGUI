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


is.windows <- function() {
    Sys.info()["sysname"] == "Windows"
}

is.mac <- function() {
    Sys.info()["sysname"] == "Darwin"
}

is.linux <- function() {
    Sys.info()["sysname"] == "Linux"
}

#' busyIndicator
#'
#' This is a function to indicate the work is in progress, it was created for the plots
#' that rendering were taking long and withprogress was not working.
#' This funcion  was adapetd from:
#' https://github.com/davesteps/homebrewR
#' and
#' https://github.com/AnalytixWare/ShinySky
#' @param text The text to show
#' @param wait The amount of time to wait before showing the busy indicator. The
#'   default is 1000 which is 1 second.
#'
#' @export
busyIndicator <- function(text = "Calculation in progress..", wait=1000) {
    tagList(
        singleton(tags$head(
            tags$link(rel="stylesheet", type="text/css",href="TCGAbiolinksGUI.css")
        )),
        div(class = "busy",
              h4(text),
              h2(HTML('<i class="fa fa-cog fa-spin"></i>'))),
              tags$script(sprintf(
                    "	setInterval(function(){
                            if ($('html').hasClass('shiny-busy')) {
                                setTimeout(function() {
                                    if ($('html').hasClass('shiny-busy')) {
                                            $('div.busy').show()
                                    }
                                }, %d)
                            } else {
                                $('div.busy').hide()
                            }
                        },100)
            ",wait)
        )
    )
}
