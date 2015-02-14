#' Create query from UI to be used with esearch - not completed
#' Create query from UI to be used with esearch - not completed
#'
#' @param List of parameters parameter = (Term,Fild,Connector)
#' @keywords internal
eGeo <- function(...) {
  input_list <- list(...)
  output_list <- lapply (X = input_list, function(x) {
    if(length(x) > 2){
      paste( x[[3]], paste( paste (x[[1]], x[[2]], sep = "["), NULL, sep = "]"))
    }
    else{
      paste( paste( NULL, x[[1]], sep = "("), NULL, sep = ")")
    }
  })
  output_list <- paste(output_list,collapse=" ")
  return(output_list)
}

#' @title Download data from GEO database
#' @description Download data from GEO database using reutils library
#' @param iQuery - GEO query - build it with http://www.ncbi.nlm.nih.gov/gds/advanced
#' @param iOut - path to save files
#' @examples
#' \dontrun{
#'    geoDownloader ("'(h1[All Fields] AND RRBS[All Fields])
#'                    AND roadmap epigenomics[Project] AND "gsm"[Filter]',
#                     "path_to_download_folder"
#'                    )
#' }
#' @seealso \url{http://www.ncbi.nlm.nih.gov/geo/info/download.html}
#' @seealso \url{http://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html}
#' @seealso \url{http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch}
#' @name geoDownloader
geoDownloader <- function(iQuery, iOut)
  {
  # Constant parameters
  roadmapURL  <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term="
  optionRet   <- "retmax=5000"
  optionHist  <- "usehistory=y"

  pmids   <- esearch (iQuery, "gds", retmax = 10000)
  if(pmids@.xData$no_errors()){
    info <- content (esummary(pmids), "parsed")
    nbFiles <<- length(info)
    ftps    <- lapply  (info, function(x)  x$FTPLink)
    dir.create (iOut, showWarnings = FALSE)

    nbFilesDownloaded  <<- 0
    for(ftp in ftps){
      dirName <- tail (unlist (strsplit (ftp,"/")), n = 1)
      outPath <- paste (iOut, dirName , sep = "/")
      dir.create (outPath, showWarnings = FALSE)
      filePath <- unlist (strsplit (getURL (paste0 (ftp, "suppl/"))," "))
      fileName <- strsplit (tail (filePath, n = 1),"\n")[[1]]
      #download.file (paste0 (ftp,"suppl/",fileName), paste0 (dirName, fileName))
    }
    print("Downloaded")
  } else {
    print("No file found")
  }
}

#' Calls UI interface
#' Calls UI interface
#' @keywords internal
biOMICsApp <- function() {
  shinyApp(  server = biOMICsServer, ui = biOMICsUI)
}
