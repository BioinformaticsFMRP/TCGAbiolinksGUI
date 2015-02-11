#' @title Download data from ENCODDE project
#' @description Downloads data from ENCODE project
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
#' @name encodeDownloader
library(RCurl)
library(rjson)

formatTarget <- function(iTarget){
  if(!is.null(iTarget)){
    return (paste(paste ("target.investigated_as", gsub(" ","%20",iTarget),sep = "="), collapse = "&"))
  }
  return(NULL)
}

formatSample <- function(iSample){
  if(!is.null(iSample)){
    return (paste(paste ("replicates.library.biosample.biosample_type", gsub(" ", "%20", iSample), sep = "="), collapse = "&"))
  }
  return (NULL)
}

formatFType <- function(iFType){
  if(!is.null(iFType)){
    return (paste(paste ("files.file_format", iFType, sep = "="), collapse = "&"))
  }
  return (NULL)
} 

formatAssay <- function(iAssay){
  if(!is.null(iAssay)){
    return (paste(paste ("assay_term_name", iAssay, sep = "="), collapse = "&"))
  }
  return (NULL)
}

formatAssembly <- function(iAssembly){
  if(!is.null(iAssembly)){
    return (paste(paste ("assembly", iAssembly, sep = "="), collapse = "&"))
  }
  return(NULL)
}

formatSearch <- function(iSearch){
  if(!is.null(iSearch)){
    return (paste ("searchTerm",  gsub(" ","+",iSearch), sep="=")) 
  }
  return (NULL)
}

formatType <- function(iType){
  if(!is.null(iType)){
    return (paste ("type", iType, sep="="))
  }
  return (NULL)
}

encodeDownloader <- function(iSearch, iType, iTarget, iSample, iFType, iAssay, iAssembly, iOut) {
  
  # Constant parameters
  encodePath  = "https://www.encodeproject.org"
  searchPath  = "/search/"
  format      = "format=json"
  frame       = "frame=embedded"
  
  searchTerm = formatSearch   ( iSearch  ) 
  type       = formatType     ( iType    )
  target     = formatTarget   ( iTarget  )
  biosample  = formatSample   ( iSample  )
  fType      = formatFType    ( iFType   )  
  assay      = formatAssay    ( iAssay   )  
  assembly   = formatAssembly ( iAssembly)

  URLoptions <- paste (searchTerm, format, frame, target, type, biosample, fType, assay, assembly, sep="&")
  URL        <- paste (encodePath, paste( searchPath, URLoptions, sep="?"), sep="")
  
  data       <- fromJSON (
    getURL (URL, dirlistonly = TRUE), 
    unexpected.escape = "keep"
  )
  
  # Output folder
  dir.create (iOut, showWarnings = FALSE)
  
  # Downloads all files from the search
  for (j in 1:length (data$'@graph')){ 
    nbFiles <- length (data$'@graph'[[j]]$files) 
    print(paste("Downloading experiment ", j, sep = ""))
    
    for (i in 1:nbFiles){
      filePath <- data$'@graph'[[j]]$files[[i]]$href
      curl = getCurlHandle()
      getURL (paste (encodePath, filePath, sep = ""), curl = curl)
      dn <- getCurlInfo(curl)$redirect.url # handling server redicrection
      
      # preparing to download - create project folder
      fileOut <- paste( iOut, unlist( strsplit( data$'@graph'[[j]]$'@id', "/"))[3], sep = "/")
      dir.create (fileOut, showWarnings = FALSE)
      fileOut <- paste( fileOut, unlist( strsplit( filePath, "/"))[5], sep = "/")
      
      # Did the user select a file format  
      if(is.null(iFType) || data$'@graph'[[j]]$files[[i]]$file_format %in% iFType){
        download.file( dn, fileOut, quiet = TRUE, method = "curl")
        print(paste("Downloaded file: ", fileOut, sep = ""))
      }

      rm(curl) 
    }
  }
}
