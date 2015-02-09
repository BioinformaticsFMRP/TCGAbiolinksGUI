#' @title Download data from ENCODDE project
#' @description
#' @param x A number
#' @return TRUE if \code{x} is prime, FALSE otherwise
#' @examples
#' prime(4)
#' prime(13)
#' \dontrun{
#'    prime('a')
#' }
#' @seealso \url{https://www.encodeproject.org/search/}
#' @seealso \url{https://www.encodeproject.org/help/rest-api/}
#' 


library(RCurl)
url <- "ftp://hgdownload.cse.ucsc.edu/goldenPath/"
filenames <- getURL(url,ftp.use.epsv = FALSE,dirlistonly = TRUE) 
filenames = paste(url, strsplit(filenames, "\n")[[1]], sep = "")


library(RJSONIO)
library(rjson)

# FIXED
encodeProject = "https://www.encodeproject.org"
encodeSearch  = "/search/"
returnParam   = "format=json"
frameParam    = "frame=embedded"

# USER SELECTION
searchParams = "searchTerm=bone+chip"
type         = "type=experiment"
target       = "target.investigated_as=transcription%20factor"

URLoptions <- paste (searchParams, returnParam, frameParam, target, type, sep="&")
URL        <- paste (encodeProject, paste( encodeSearch, URLoptions, sep="?"), sep="")

filenames  <- getURL (URL, dirlistonly = TRUE) 
data       <- fromJSON (filenames, unexpected.escape = "keep")

# Downloading folder
dir.create("encodeDownload")

  
# Downloads all files from the search
for (j in 1:length (data$'@graph')){ 
  for (i in 1:length (data$'@graph'[[j]]$files)){
    aux <- data$'@graph'[[j]]$files[[i]]$href
    dn <- paste (encodeProject, aux, sep = "")
    curl = getCurlHandle()
    getURL (dn, curl = curl)
    dn <- getCurlInfo(curl)$redirect.url
    rm(curl) # release the curl!
    download.file( dn, paste("encodeDownload",unlist(strsplit(aux,"/"))[5],sep="/"), method="curl")
  }
}
