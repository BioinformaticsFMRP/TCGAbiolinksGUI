#' @title biOMICs.download
#' @description Download data previously selected using the biomics.search function
#' @param lines biomics.search output
#' @param enc.file.type Extension to be downloaed from encode database
#' @param rmap.file.type Extension to be downloaed from roadmap database
#' @seealso biOMICs.search
#' @export
biOMICs.download <- function(lines,
                             enc.file.type = NULL,
                             rmap.file.type=NULL
){

  encode.lines <- subset(lines, database == 'encode')
  rmap.lines <- subset(lines, database == 'roadmap')
  tcga.lines <- subset(lines, database == 'tcga')

  #--------------   download Encode
  if(dim(encode.lines)[1] >0){
    ENCODEDownload(encode.lines,enc.file.type)
  }

  #--------------   download ROADMAP
  if(dim(rmap.lines)[1] > 0){
    ROADMAPDownload(rmap.lines, rmap.file.type)
  }
  #---------------- download TCGA
  # TODO: add filters, folder to save
  if(dim(tcga.lines)[1] > 0){
    tcga.lines <- tcga.db[tcga.db$name == tcga.lines$ID,]
    TCGADownload(tcga.lines,path="TCGA")
  }
}

ENCODEDownload <- function(lines,type = NULL){
  encode.url <- "https://www.encodeproject.org/"
  json <- "/?format=JSON"
  for (i in 1:dim(lines)[1]){
    id <-  lines$ID[i]
    url <- paste0(encode.url, 'experiments/', id, json)

    #print(url)
    #get list of files
    item <- rjson::fromJSON (
      RCurl::getURL(url, dirlistonly = TRUE,
                    .opts = list(ssl.verifypeer = FALSE))
    )[['files']]

    files <- sapply(item,function(x){x$href})

    if(!is.null(type)){
      idx <- unlist(lapply(type,function(x){grep(x,files)}))
      files <- files[idx]
    }

    #download files
    for (j in seq_along(item)){
      file <- item[[j]]$href
      link <- paste0(encode.url, file)

      fileout <- unlist(strsplit(link, "/"))
      fileout <- fileout[length(fileout)]

      downloader::download(link, fileout)
    }
  }

}

#' @import RCurl
ROADMAPDownload <- function (lines,
                             type=NULL,
                             path="."
){

  dir.create(path,showWarnings = F)

  for (i in 1:dim(lines)[1]){
    id <-  lines$ID[i]
    url <- roadmap.db[roadmap.db$X..GEO.Accession == id,]$GEO.FTP
    if(length(url) == 0){next}
    #get list of files in roadmap FTP
    filenames <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
    filenames <- strsplit(filenames, "\r\n")
    filenames <- unlist(filenames)

    if(!is.null(type)){
      idx <- unlist(lapply(type,function(x){grep(x,filenames)}))
      filenames <- filenames[idx]
    }
    files <- paste0(path,"/",basename(filenames[j]))
    #download files
    print(files)
    for(j in seq_along(files)){
      print(url)
      print(files)

      #downloader::download(paste0(url,filenames[j]), files[j])
    }
  }
}

#' @title Download roadmap data
#' @description Download roadmap data previously selected using the roadmap.search function
#' @param lines roadmap.search output
#' @param path Folder to save the file
#' @export
#' @import RCurl
get.roadmap <- function (lines,
                         type=NULL,
                         path="."
){

  dir.create(path, showWarnings = F, recursive =T)

  for (i in 1:dim(lines)[1]){
    url <- lines[i,]$GEO.FTP
    if(length(url) == 0){next}
    filenames <- getFileNames(url)
    if(!is.null(type)){
      idx <- unlist(lapply(type,function(x){grep(x,filenames)}))
      filenames <- filenames[idx]
    }
    files <- paste0(url,filenames)

    #download files
    for(j in seq_along(files)){
      aux <- paste0(path,"/",basename(files[j]))
      if(!file.exists(aux)){
        downloader::download(files[j],aux)
      }
    }
  }
}
#' @title TCGA Download
#' @description Download data previously selected using the TCGAQuery function
#' @param data TCGAQUery output
#' @param path location of the final data saving
#' @seealso TCGAQuery
#' @examples
#' \dontrun{
#'          TCGADownload(data,"folder")
#' }
#' @export
#' @import downloader
TCGADownload <- function(data=NULL,path=".")
{
  dir.create(path,showWarnings = F)
  root <- "https://tcga-data.nci.nih.gov"
  if(!("file" %in% colnames(data))){
    message("Downloading folders")
    for(i in 1:nrow(data)){
      file <- paste0(path,"/",basename(data[i,"deployLocation"]))
      message(paste0("Downloading:", basename(data[i,"deployLocation"])))
      if(!file.exists(file)){
        downloader::download(paste0(root,data[i,"deployLocation"]),file)
        untar(file, exdir = path)
      }
    }
  }
  else{
    message("Downloading files")
    for(i in 1:nrow(data)){
      file <- paste0(path,"/",basename(data[i,"file"]))
      message(paste0("Downloading:", basename(data[i,"file"])))
      if(!file.exists(file)){
        downloader::download(paste0(root,gsub(".tar.gz","",data[i,"deployLocation"]),"/",data[i,"file"]),file)
      }
    }
  }
}
