#' term=RRBS%5BAll%20Fields%5D%20AND%20H1%5BAll%20Fields%5D%20AND%20%22gsm%22%5BFilter%5D%20&cmd=DetailsSearch
#' RRBS[All Fields] AND H1[All Fields] AND roadmap epigenomics[Project]
# http://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html
# X[[3]] = Connector (AND,OR, NOT)
# X[[2]] = Field [Project]
# X[[1]] = Term "Rattus norvegicus"
searchGeo <- function(...) {
  input_list <- list(...)
  output_list <- lapply(X=input_list, function(x) {
    if(length(x) > 1){
      paste(paste(x[[1]],gsub(" ","%20",x[[2]]),sep="["),NULL,sep="]")
    }
    else{
      paste(paste(NULL,x[[1]],sep="("),NULL,sep=")")
    }
  })
  output_list <- paste(output_list,collapse="")
  return(output_list)
}

eGeo <- function(...) {
  input_list <- list(...)
  output_list <- lapply(X=input_list, function(x) {
    if(length(x) > 1){
      paste(x[[3]],paste(paste(x[[1]],x[[2]],sep="["),NULL,sep="]"))
    }
    else{
      paste(paste(NULL,x[[1]],sep="("),NULL,sep=")")
    }
  })
  output_list <- paste(output_list,collapse="")
  return(output_list)
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
#' 


# GSE[ETYP]+AND+"published+last+3+months"[Filter]

geoDownloader <- function(iSearch, iType, iTarget,
                             iSample, iFType, iAssay,
                             iAssembly, iOut)
  {
  
  # Constant parameters
  roadmapURL  <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term="
  optionRet    <- "retmax=5000"
  optionHist   <- "usehistory=y"
  
  search <- searchGeo(c("h1 RRBS",NULL,NULL),c("roadmap epigenomics","Project","AND"),c("gsm","Filter","AND"))

  query <- '(h1[All Fields] AND RRBS[All Fields]) AND roadmap epigenomics[Project] AND "gsm"[Filter]'
  pmids <- esearch(query, "gds")
  content <- content(esummary(pmids), "parsed")
  ftps <- lapply(content, function(x)  x$FTPLink) 
  dir.create ("downloadRM", showWarnings = FALSE)

  for(ftp in ftps){
    dirName <- unlist( strsplit(ftp,"/"))
    dirName <- paste("downloadRM",dirName[length(dirName)],sep="/")
    dir.create (dirName, showWarnings = FALSE)
    file <- unlist( strsplit(getURL(paste(ftp,"suppl/",sep=""))," "))
    file <- file[length(file)]                          
    fileName <- strsplit(file,"\n")[[1]]
    download.file(paste(ftp,"suppl/",fileName,sep=""), paste(dirName,fileName,sep="/"))
  }
  print("Downloaded")  
 
}

#' Calls UI interface
#' Calls UI interface
#' @keywords internal
roadmapApp <- function() {
  shinyApp(  server = roadmapServer, ui = roadmapUI)
}
