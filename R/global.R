#' @title biOmics interface
#' @description Calls UI interface
#' @examples
#'    biOmicsApp()
#' @seealso \url{https://www.encodeproject.org/search/}
#' @seealso \url{https://www.encodeproject.org/help/rest-api/}
#' @name biOmicsApp
#' @import shiny shinyFiles shinydashboard SummarizedExperiment
#' TCGAbiolinks UpSetR ggplot2 shinyBS stringr ggrepel pathview ELMER grid
#' @export
#' @return Open a connection to shiny
biOmicsApp <- function() {
    shiny::runApp(system.file("app", package = "biOmics"),launch.browser=TRUE)
}


#' @title TCGAbiolinks GUI
#' @description Calls UI interface
#' @examples
#'    TCGAbiolinksGUI()
#' @name biOmicsApp
#' @import shiny shinyFiles shinydashboard SummarizedExperiment
#' TCGAbiolinks UpSetR ggplot2 shinyBS stringr ggrepel pathview ELMER grid
#' @export
#' @return Open a connection to shiny
TCGAbiolinksGUI <- function() {
    shiny::runApp(system.file("app", package = "biOmics"),launch.browser=TRUE)
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

#' @export
get.obj <- function(obj){
    return(get(obj))
}
