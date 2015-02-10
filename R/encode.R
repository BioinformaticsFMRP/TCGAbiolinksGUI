#' @title Download data from ENCODDE project
#' @description
#' @param searchTerm - bone chip
#' @param type - "experiment" "publication" "software" ""antibody_lot" "web"
#' @param target - "transcription factor" "RNA binding protein" "tag" "histone modification"
#' @param output - path to save files
#' @examples
#' encodeDownload("bone chip",
#'                "experiment",
#'                "transcription factor", 
#'                "encodeDownload"
#'                )
#' @seealso \url{https://www.encodeproject.org/search/}
#' @seealso \url{https://www.encodeproject.org/help/rest-api/}
#' 
library(RCurl)
library(rjson)

# USER SELECTION - needs improvement
search  = "bone chip"
type    = "experiment"
target  = "control" # can be more than one
output  = "encodeDownload"

iSearch = search
iType= type
iTarget = target
iOut = output


encodeDownload <- function(iSearch, iType, iTarget, iOut) {

  # Constant parameters
  encodePath  = "https://www.encodeproject.org"
  searchPath  = "/search/"
  format      = "format=json"
  frame       = "frame=embedded"
  
  searchTerm  = paste ("searchTerm",  gsub(" ","+",iSearch), sep="=")
  type        = paste ("type", iType, sep="=")
  target      = paste ("target.investigated_as", gsub(" ","%20",iTarget),sep="=")
  
  URLoptions <- paste (searchTerm, format, frame, target, type, sep="&")
  URL        <- paste (encodePath, paste( searchPath, URLoptions, sep="?"), sep="")
  
  data       <- fromJSON (
                  getURL (URL, dirlistonly = TRUE), 
                  unexpected.escape = "keep"
                )

  # Output folder
  dir.create(iOut)
  
  # Downloads all files from the search
  for (j in 1:length (data$'@graph')){ 
      nbFiles <- length (data$'@graph'[[j]]$files) 
      print(paste("Experiment ",j,": Downloading ", nbFiles, " files", sep = ""))
    for (i in 1:nbFiles){
      filePath <- data$'@graph'[[j]]$files[[i]]$href
      curl = getCurlHandle()
      getURL (paste (encodePath, filePath, sep = ""), curl = curl)
      dn <- getCurlInfo(curl)$redirect.url # handling server redicrection
      
      # preparing to download
      fileOut <- paste( iOut, unlist( strsplit( filePath, "/"))[4], sep = "/")
      dir.create(fileOut)
      fileOut <- paste( iOut, unlist( strsplit( filePath, "/"))[5], sep = "/")
      #download.file( dn, fileOut, method = "curl")
      download.file( dn, fileOut, quiet = TRUE, method = "auto")
      rm(curl) 
    }
  }
}