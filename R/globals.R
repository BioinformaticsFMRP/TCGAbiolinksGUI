# @title Reactive global variables
# @description Global variables used to update shiny
# @keywords internal
.result <- reactiveValues(g_nbFiles = 0, g_downloaded = 0, df = NULL, downloading = 1)

#' @title biOMICs interface
#' @description Calls UI interface
#' @examples
#' \dontrun{
#'    biOMICsApp()
#' }
#' @seealso \url{https://www.encodeproject.org/search/}
#' @seealso \url{https://www.encodeproject.org/help/rest-api/}
#' @name biOMICsApp
#' @keywords internal
biOMICsApp <- function() {
  shiny::runApp (system.file ('app', package = 'biOMICs'))
}
