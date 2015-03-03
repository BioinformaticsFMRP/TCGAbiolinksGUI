# @title Reactive global variables
# @description Global variables used to update shiny
# @keywords internal
result <- reactiveValues(g_nbFiles = 0, g_downloaded = 0, df = NULL, downloading = 1)

.onAttach <- function (libname, pkgname){
  updateRoadmap()
  load(file = system.file("extdata/roadmap.RData", package="biOMICs"), .GlobalEnv)
}

updateRoadmap <- function (){
  #get raodmap data
  roadmapSummaryFile <- system.file("extdata/roadmap.RData", package="biOMICs")
  dataDirecotry <- system.file("extdata", package="biOMICs")

  # If file exists and it was downlaoded 7 day ago ignore
  # otherwise redownload it
  if (file.exists(roadmapSummaryFile)){
    finf <- file.info(roadmapSummaryFile, extra_cols = FALSE)
    if (nrow(finf[difftime(Sys.time(), finf[,"mtime"], units = "days") > 7 , 1:4]) == 2){
      print("file needs update")
      downloader::download("http://www.ncbi.nlm.nih.gov/geo/roadmap/epigenomics/?view=samples&sort=acc&mode=csv",
                           paste0(dataDirecotry,"/roadmap.csv"))
      roadmap <- NULL
      roadmap$summary        <- read.csv(roadmapSummaryFile)
      roadmap$experimentList <- levels(unique(roadmap$summary$Experiment))
      roadmap$samplesList    <- levels(unique(roadmap$summary$Sample.Name))
      save(roadmap,file = paste0(dataDirecotry,"/roadmap.RData"))
      file.remove (paste0(dataDirecotry,"/roadmap.csv"))
    }
  }


}

#' @title biOMICs interface
#' @description Calls UI interface
#' @examples
#' \dontrun{
#'    biOMICsApp()
#' }
#' @seealso \url{https://www.encodeproject.org/search/}
#' @seealso \url{https://www.encodeproject.org/help/rest-api/}
#' @name biOMICsApp
#' @import shiny shinydashboard
#' @export
biOMICsApp <- function() {
  shiny::runApp (system.file ('app', package = 'biOMICs'))
}
