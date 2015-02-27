#' Create query from UI to be used with esearch - not completed
#' Create query from UI to be used with esearch - not completed
#' @param List of parameters parameter = (Term,Field,Connector)
#'        First parameter should contain the terms,
#'        Field NULL which means all fields and connector NULL
#'        For roadmap one parameter should be c ("roadmap epigenomics","Project","AND"),
#'        For samples: c ("gsm","Filter","AND") )
#' @examples
#' \dontrun{
#'  eGeo(
#'    c ("h1 cell RRBS",NULL,NULL),
#'    c ("roadmap epigenomics","Project","AND"),
#'    c ("gsm","Filter","AND") )
#' }
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
#'              Input should be a GEO query, you can you this site
#'              to create the query  http://www.ncbi.nlm.nih.gov/gds/advanced
#' @param iQuery - GEO query - build it with http://www.ncbi.nlm.nih.gov/gds/advanced
#' @param iOut - path to save files
#' @examples
#' \dontrun{
#'    geoDownloader ("'(h1[All Fields] AND RRBS[All Fields])
#'                    AND roadmap epigenomics[Project] AND "gsm"[Filter]',
#                     "path_to_download_folder"
#'                    )
#'  iQuery <- "((GSM956006) AND gsm[Filter]) AND roadmap epigenomics[Project]"
#'  iQuery <- "((h1 cell RRBS) AND gsm[Filter]) AND roadmap epigenomics[Project]"
#'  iQuery <- "((h1 cell ) AND gsm[Filter]) AND roadmap epigenomics[Project]"
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

  # search on geo databse(gds) for iQuery
  pmids <- reutils::esearch (iQuery, "gds", retmax = 10000)

  if (pmids@.xData$no_errors()) {
    nbFiles  <- as.numeric (XML::xmlToList (pmids@.xData$get_content())$Count)

    if (nbFiles  > 0 ){ # Found more than one result?
      .mkdir(iOut)

      suppressWarnings({
        info <- reutils::content (reutils::esummary(pmids), "parsed") # get files info
      })

      if (nbFiles >1) ftps <- lapply  (info, function(x)  x$FTPLink)
      else ftps <- info$FTPLink

      #ftpsSRA <- lapply  (info$ ExtRelations.ExtRelation.TargetFTPLink, function(x)  x)

      for(ftp in ftps){
        if(!is.na(ftp)){

          # create folder for experiment
          dirName <- tail (unlist (strsplit (ftp,"/")), n = 1) # get experiment Name
          outPath <- paste (iOut, dirName , sep = "/")
          .mkdir (outPath)

          # get list of files names in the ftp
          ftp <- paste0 (ftp, "suppl/")
          fileName <- .getFileNames(ftp)
          ret <- .downloadFromGEO(ftp,fileName,outPath)
        }
      }

      .prepareInfoTable(ret$link,iOut)
      .uncompress(ret$compresssedFiles)

    }
  }
}
.downloadFromGEO <- function(iFTP, iFileName,iOut){
  link <- list()
  compresssedFiles <- c()
  lapply (iFileName,
          function(x){
            ftpLink <- paste0 (iFTP, x)
            file    <- paste0 (iOut, "/", x)
            if (!file.exists(file)){
              downloader::download(ftpLink, file)
              compresssedFiles <- c(compresssedFiles, file)
            }
            link <- c(link, ftpLink)
          }
  )
  ret <- NULL
  ret$link <- link
  ret$compresssedFiles <- compresssedFiles

}

geoDownloaderLinks <- function(iLinks, iOut)
{

  nbFiles  <- length(iLinks)
  if (nbFiles  > 0 ){ # Found more than one result?
    .mkdir(iOut)

    for(ftp in iLinks){
      if(!is.na(ftp)){
        # create folder for experiment
        dirName <- tail (unlist (strsplit (ftp,"/")), n = 1) # get experiment Name
        outPath <- paste (iOut, dirName , sep = "/")
        .mkdir (outPath)
        # get list of files names in the ftp
        fileName <- .getFileNames(ftp)
        ret <- .downloadFromGEO(ftp,fileName,outPath)
      }
    }
    .prepareInfoTable(ret$link,iOut)
    .uncompress(ret$compresssedFiles)
  }
}

# Uncompress a list of files - it does not remove the compressed file
# @keywords internal
.uncompress <- function(iFiles){
  if (!is.null(iFiles)){
    lapply  (iFiles,
             function(x){
               if(file.exists(x))
                 R.utils::gunzip(x, remove = FALSE)
             }
    )
  } else {
    return(NULL)
  }
}

# Show ftp links downloaded
# @keywords internal
.prepareInfoTable <- function(iLink,iOut){
  if(length(iLink) > 0 ){
    df <- do.call (rbind.data.frame, iLink)
    .result$g_nbFiles <- nrow(df)
    names(df) <- paste0 ("Files downloaded into:", getwd(),"/",iOut)
    .result$df <- df
  }
}

# Get all files in the ftp folder
# @keywords internal
.getFileNames <- function(ftp){
  print(ftp)
  filePath <- unlist (
    strsplit(
      RCurl::getURL(
        ftp,
        ftp.use.epsv = FALSE,
        dirlistonly = TRUE
      ),
      " "
    )
  )
  # remove \r\n from name
  fileName <- as.list(strsplit (tail (filePath, n = 1),"\r*\n")[[1]])
}

# Create directory
# @keywords internal
.mkdir <- function (iOut){
  if(!file.exists(iOut))
    dir.create (iOut, showWarnings = TRUE)  # create directory to save files
}
