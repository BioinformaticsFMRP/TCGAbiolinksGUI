#' Concatenate Target term for search
#' Concatenate Target term for search
#'
#' @param list to concatenate
#' @keywords internal
formatTarget <- function(iTarget){
  if( !is.null(iTarget)) {
    return (paste(
                  "target.investigated_as",
                  gsub(" ","%20",iTarget),
                  sep = "=",
                  collapse = "&"
                  )
            )
  }
  return(NULL)
}

#' Concatenate Sample term for search
#' Concatenate Sample term for search
#'
#' @param list to concatenate
#' @keywords internal
formatSample <- function(iSample){
  if (!is.null(iSample)) {
    return (paste(
                  "replicates.library.biosample.biosample_type",
                  gsub(" ", "%20", iSample),
                  sep = "=",
                  collapse = "&")
            )
  }
  return (NULL)
}

#' Concatenate Assay term for search
#' Concatenate Assay term for search
#'
#' @param list to concatenate
#' @keywords internal
formatFType <- function(iFType){
  if (!is.null(iFType)) {
    return (paste(
                  "files.file_format",
                  iFType,
                  sep = "=",
                  collapse = "&")
            )
  }
  return (NULL)
} 

#' Concatenate Assay term for search
#' Concatenate Assay term for search
#'
#' @param list to concatenate
#' @keywords internal
formatAssay <- function(iAssay){
  if( !is.null(iAssay)) {
    return (paste(
                  "assay_term_name",
                  iAssay,
                  sep = "=",
                collapse = "&")
            )
  }
  return (NULL)
}

#' Concatenate Assembly term for search
#' Concatenate Assembly term for search
#'
#' @param list to concatenate
#' @keywords internal
formatAssembly <- function(iAssembly){
  if(!is.null(iAssembly)){
    return (paste(
                  "assembly",
                  iAssembly,
                  sep = "=",
                collapse = "&")
            )
  }
  return(NULL)
}


#' Concatenate search term for search
#' Concatenate search term for search
#'
#' @param list to concatenate
#' @keywords internal
formatSearch <- function(iSearch){
  if(!is.null(iSearch))
    return (paste ("searchTerm",  gsub(" ","+",iSearch), sep="="))
  return (NULL)
}


#' Concatenate type term for search
#' Concatenate type term for search
#'
#' @param list to concatenate
#' @keywords internal
formatType <- function(iType){
  if(!is.null(iType))
    return (paste ("type", iType, sep="="))
  return (NULL)
}

#' @title Download data from ENCODDE project
#' @description Downloads data from ENCODE project
#' @param iSearch - bone chip
#' @param iType -   "experiment" "publication" "software" ""antibody_lot" "web"
#' @param iTarget - "transcription factor" "RNA binding protein" 
#'                  "tag" "histone modification"
#' @param iSample - "tissue" "primary cell"
#' @param iFType  - "bam" "bigWig" "bed_broadPeak" "broadPeak" "fastq"
#' @param iAssay - "ChIP-seq" "RIP-chip" "Repli-chip" 
#' @param iAssembly - "hg19" "mm9"
#' @param iOut - path to save files
#' @examples
#' \dontrun{
#'    encodeDownloader("bone chip",
#'                    "experiment",
#'                    "transcription factor", 
#'                    "primary cell",
#'                    "bigWig",
#'                    "ChIP-seq",
#'                    "mm9",
#'                    "encodeDownload"
#'                    )
#' }
#' @seealso \url{https://www.encodeproject.org/search/}
#' @seealso \url{https://www.encodeproject.org/help/rest-api/}
#' @name encodeDownloader
encodeDownloader <- function(iSearch, iType, iTarget,
                             iSample, iFType, iAssay,
                             iAssembly, iOut)
  {
  
  # Constant parameters
  encodePath  <- "https://www.encodeproject.org"
  searchPath  <- "/search/"
  format      <- "format=json"
  frame       <- "frame=embedded"
  
  searchTerm <- formatSearch   ( iSearch  )
  type       <- formatType     ( iType    )
  target     <- formatTarget   ( iTarget  )
  biosample  <- formatSample   ( iSample  )
  fType      <- formatFType    ( iFType   )
  assay      <- formatAssay    ( iAssay   )
  assembly   <- formatAssembly ( iAssembly)
  
  URLoptions <- paste (searchTerm, format, frame, target,
                       type, biosample, fType, assay, assembly, 
                       sep="&")
  URL        <- paste (encodePath, 
                       paste( searchPath, URLoptions, sep="?"),
                       sep="")
  
  data       <- rjson::fromJSON (
    RCurl::getURL (URL, 
                   dirlistonly = TRUE,
                   .opts = list(ssl.verifypeer = FALSE)
                   ),
    unexpected.escape = "keep"
  )
  
  if (data$total > 0){
    # Output folder
    dir.create (iOut, showWarnings = FALSE)
  
    # Downloads all files from the search
    for (j in 1:length (data$'@graph')){
      nbFiles <- length (data$'@graph'[[j]]$files)
      print(paste("Downloading experiment ", j, sep = ""))
    
      for (i in 1:nbFiles){
        filePath <- data$'@graph'[[j]]$files[[i]]$href
        dn <- paste (encodePath, filePath, sep = "")
        fileOut <- paste( iOut, 
                        unlist( strsplit( data$'@graph'[[j]]$'@id', "/"))[3],
                        sep = "/")
        dir.create (fileOut, showWarnings = FALSE)
        fileOut <- paste( fileOut, 
                          unlist( strsplit( filePath, "/"))[5], 
                          sep = "/")
      
        # Did the user select a file format  
        if(is.null(iFType) || 
             data$'@graph'[[j]]$files[[i]]$file_format %in% iFType
          ){
          download( dn, fileOut, quiet = TRUE)
          print(paste("Downloaded file: ", fileOut, sep = ""))
        }
      
      }
    }
    print("Downloaded")  
  } else { 
      print(data$notification)
  }
}

#' Calls UI interface
#' Calls UI interface
#' @keywords internal
encodeUi <- function() {
  shinyApp(  server = shinyServer, ui = shinyUI)
}