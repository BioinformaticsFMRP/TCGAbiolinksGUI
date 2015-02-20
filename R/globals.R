# @title Reactive global variables
# @description Global variables used to update shiny
# @keywords internal
.result <- reactiveValues(g_nbFiles = 0, g_downloaded = 0, df = NULL, downloading = 1)
