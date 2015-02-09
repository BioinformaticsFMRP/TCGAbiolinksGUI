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
searchTerm  = "bone chip"
type        = "experiment"
target      = "transcription factor" # can be more than one
output      = "encodeDownload"

encodeDownload <- function(search, type, target, output) {

  # Constant parameters
  encodePath  = "https://www.encodeproject.org"
  searchPath  = "/search/"
  format      = "format=json"
  frame       = "frame=embedded"
  
  searchTerm  = paste ("searchTerm",  gsub(" ","+",search), sep="=")
  type        = paste ("type", type, sep="=")
  target      = paste ("target.investigated_as", target,sep="=")
  
  URLoptions <- paste (searchTerm, format, frame, target, type, sep="&")
  URL        <- paste (encodePath, paste( searchPath, URLoptions, sep="?"), sep="")
  
  data       <- fromJSON (
                  getURL (URL, dirlistonly = TRUE), 
                  unexpected.escape = "keep"
                )

  # Output folder
  dir.create(output)
  
  # Downloads all files from the search
  for (j in 1:length (data$'@graph')){ 
    for (i in 1:length (data$'@graph'[[j]]$files)){
      curl = getCurlHandle()
      dn <- paste (encodeProject, data$'@graph'[[j]]$files[[i]]$href, sep = "")
      getURL (dn, curl = curl)
      dn <- getCurlInfo(curl)$redirect.url
      rm(curl) # release the curl!
      download.file( dn, paste("encodeDownload",unlist(strsplit(aux,"/"))[5],sep="/"), method="curl")
      rm(dn)
    }
  }
}