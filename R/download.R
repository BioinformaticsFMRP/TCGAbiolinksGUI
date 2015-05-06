#' @export
biOMICs.download <- function(lines){

  encode.url <- "https://www.encodeproject.org/"
  tcga.url <- "https://tcga-data.nci.nih.gov"
  json <- "/?format=JSON"

  encode.lines <- subset(lines, database == 'encode')
  rmap.lines <- subset(lines, database == 'roadmap')
  tcga.lines <- subset(lines, database == 'tcga')

  #--------------   download Encode
  if(dim(encode.lines)[1] >0){
    ENCODEDownload(encode.lines)
  }

  #--------------   download ROADMAP
  if(dim(rmap.lines)[1] > 0){
    ROADMAPDownload(rmap.lines)
  }
  #---------------- download TCGA
  # TODO: add filters, folder to save
  if(dim(tcga.lines)[1] > 0){
    tcga.lines <- tcga.db[tcga.db$name == tcga.lines$ID,]
    TCGADownload(tcga.lines,path="TCGA")
  }
}

ENCODEDownload <- function(lines){
  for (i in 1:dim(encode.lines)[1]){
    id <-  encode.lines$ID[i]
    url <- paste0(encode.url, 'experiments/', id, json)

    print(url)
    #get list of files
    item <- rjson::fromJSON (
      RCurl::getURL(url, dirlistonly = TRUE,
                    .opts = list(ssl.verifypeer = FALSE))
    )[['files']]

    #download files
    for (j in 1:length(item)){
      file <- item[[j]]$href
      link <- paste0(encode.url, file)

      fileout <- unlist(strsplit(link, "/"))
      fileout <- fileout[length(fileout)]

      print(link)
      downloader::download(link, fileout)
    }
  }

}
ROADMAPDownload <- function (lines){

  for (i in 1:dim(lines)[1]){
    id <-  roadmap$ID[i]
    url <- roadmap.db[roadmap.db$X..GEO.Accession == id,]$GEO.FTP

    #get list of files in roadmap FTP
    filenames <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
    filenames <- strsplit(filenames, "\r\n")
    filenames <- unlist(filenames)

    #download files
    for(j in 1:length(filenames)){
      downloader::download(paste0(url,filenames[j]), filenames[j])
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
