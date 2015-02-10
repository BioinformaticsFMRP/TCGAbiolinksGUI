#' @title Download data from ENCODDE project
#' @description
#' @param searchTerm - bone chip
#' @param iType -   "experiment" "publication" "software" ""antibody_lot" "web"
#' @param iTarget - "transcription factor" "RNA binding protein" "tag" "histone modification"
#' @param iSample - "tissue" "primary cell"
#' @param iFType  - "bam" "bigWig" "bed_broadPeak" "broadPeak" "fastq"
#' @param assay - "ChIP-seq" "RIP-chip" "Repli-chip" 
#' @param assembly - "hg19" "mm9"
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
search    = "bone chip"
type      = "experiment"
#target    = c("control","transcription factor") # can be more than one
target    = c("control")
biosample = c("primary cell")
#biosample = c("tissue", "primary cell")
fType     = c("bigWig")
fType     = NULL
assay     = c("ChIP-seq")
assembly  = c("mm9")
output    = "encodeDownload"

iSearch   = search
iType     = type
iTarget   = target
iOut      = output
iSample   = biosample
iFType    = fType 
iAssembly = assembly
iAssay    = assay

encodeDownload <- function(iSearch, iType, iTarget, iSample, iFType, iAssay, iAssembly, iOut) {

  # Constant parameters
  encodePath  = "https://www.encodeproject.org"
  searchPath  = "/search/"
  format      = "format=json"
  frame       = "frame=embedded"
  fType       = NULL
  target      = NULL
  biosample   = NULL
  assay       = NULL
  assembly    = NULL
  
  searchTerm  = paste ("searchTerm",  gsub(" ","+",iSearch), sep="=")
  type        = paste ("type", iType, sep="=")

  if(!is.null(iTarget)){
    target = paste(paste ("target.investigated_as", gsub(" ","%20",iTarget),sep = "="), collapse = "&")
  }
  
  if(!is.null(iSample)){
    biosample = paste(paste ("replicates.library.biosample.biosample_type", gsub(" ", "%20", iSample), sep = "="), collapse = "&")
  }
  
  if(!is.null(iFType)){
    fType = paste(paste ("files.file_format", iFType, sep = "="), collapse = "&")
  }
  
  if(!is.null(iAssay)){
    assay = paste(paste ("assay_term_name", iAssay, sep = "="), collapse = "&")
  }
  if(!is.null(iAssembly)){
    assembly = paste(paste ("assembly", iAssembly, sep = "="), collapse = "&")
  }
  
  URLoptions <- paste (searchTerm, format, frame, target, type, biosample, fType, assay, assembly, sep="&")
  URL        <- paste (encodePath, paste( searchPath, URLoptions, sep="?"), sep="")
  
  data       <- fromJSON (
                  getURL (URL, dirlistonly = TRUE), 
                  unexpected.escape = "keep"
                )

  # Output folder
  dir.create(iOut, showWarnings = FALSE)
  
  # Downloads all files from the search
  for (j in 1:length (data$'@graph')){ 
      nbFiles <- length (data$'@graph'[[j]]$files) 
      print(paste("Downloading experiment ",j, sep = ""))
    for (i in 1:nbFiles){
      filePath <- data$'@graph'[[j]]$files[[i]]$href
      curl = getCurlHandle()
      getURL (paste (encodePath, filePath, sep = ""), curl = curl)
      dn <- getCurlInfo(curl)$redirect.url # handling server redicrection
      
      # preparing to download - create project folder
      fileOut <- paste( iOut, unlist( strsplit( data$'@graph'[[j]]$'@id', "/"))[3], sep = "/")
      dir.create(fileOut, showWarnings = FALSE)
      fileOut <- paste( fileOut, unlist( strsplit( filePath, "/"))[5], sep = "/")

      # Did my user select a file format  
      if(is.null(iFType) || data$'@graph'[[j]]$files[[i]]$file_format %in% iFType){
        download.file( dn, fileOut, quiet = TRUE, method = "curl")
        print(paste("Downloaded file: ", fileOut, sep = ""))
      }
      
      rm(curl) 
    }
  }
}
