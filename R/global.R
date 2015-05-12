#' @title biOmics interface
#' @description Calls UI interface
#' @examples
#'    biOmicsApp()
#' @seealso \url{https://www.encodeproject.org/search/}
#' @seealso \url{https://www.encodeproject.org/help/rest-api/}
#' @name biOmicsApp
#' @import shiny shinydashboard tcltk
#' @export
#' @return Open a connection to shiny
biOmicsApp <- function() {
    shiny::runApp(system.file("app", package = "biOmics"))
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
