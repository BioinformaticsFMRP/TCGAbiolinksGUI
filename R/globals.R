#' Reactive global variables
#' Global variables used to update shiny
result <- reactiveValues(g_nbFiles = 0, g_downloaded = 0, df = NULL, downloading = 1)
