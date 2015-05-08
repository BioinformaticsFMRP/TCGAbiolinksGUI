#' @title biOmics.download
#' @description Download data previously selected using the biOmics.search function
#' @param lines biOmics.search output
#' @param enc.file.type Extension to be downloaed from encode database
#' @param rmap.file.type Extension to be downloaed from roadmap database
#' @seealso biOmics.search
#' @export
biOmics.download <- function(lines,
                             enc.file.type = NULL,
                             rmap.file.type=NULL
){

  encode.lines <- subset(lines, database == 'encode')
  rmap.lines <- subset(lines, database == 'roadmap')
  tcga.lines <- subset(lines, database == 'tcga')

  #--------------   download Encode
  if(dim(encode.lines)[1] >0){
    encode.lines <- subset(encode.db, encode.db$accession == encode.lines$ID)
    encode.download(encode.lines,enc.file.type)
  }

  #--------------   download ROADMAP
  if(dim(rmap.lines)[1] > 0){
    rmap.lines <- subset(roadmap.db, roadmap.db$X..GEO.Accession == rmap.lines$ID)
    roadmap.download(rmap.lines,rmap.file.type)
  }
  #---------------- download TCGA
  # TODO: add filters, folder to save
  if(dim(tcga.lines)[1] > 0){
    tcga.lines <- tcga.db[tcga.db$name == tcga.lines$ID,]
    TCGADownload(tcga.lines,path="TCGA")
  }
}

#' @title Download encode data
#' @description Download encode data previously selected using the encode.search function
#' @param lines encode.search output
#' @param path Folder to save the file
#' @param type extesion of files to be downloaded
#' @export
encode.download <- function(lines,
                            type = NULL,
                            path="."
){
  dir.create(path, showWarnings = F, recursive =T)

  encode.url <- "https://www.encodeproject.org/"
  json <- "/?format=JSON"

  for (i in 1:dim(lines)[1]){
    id <-  lines$accession[i]
    url <- paste0(encode.url, 'experiments/', id, json)

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
    for (j in seq_along(files)){
      link <- paste0(encode.url, files[j])
      fileout <- file.path(path,id,basename(link))
      downloader::download(link, fileout)
    }
  }

}

#' @title Download roadmap data
#' @description Download roadmap data previously selected using the roadmap.search function
#' @param lines roadmap.search output
#' @param path Folder to save the file
#' @param type extesion of files to be downloaded
#' @export
#' @import RCurl
roadmap.download <- function (lines,
                         type=NULL,
                         path="."
){

  dir.create(path, showWarnings = F, recursive =T)
  error <- c()
  for (i in 1:dim(lines)[1]){
    url <- lines[i,]$GEO.FTP
    if(nchar(url) == 0){
      error <- c(error, lines[i,]$X..GEO.Accession)
      next
    }
    filenames <- getFileNames(url)
    if(!is.null(type)){
      idx <- unlist(lapply(type,function(x){grep(x,filenames)}))
      filenames <- filenames[idx]
    }
    files <- paste0(url,filenames)
    name.sample <- str_replace_all(lines[i,]$Sample.Name, "[^[:alnum:]]", "_")
    name.experiment <- str_replace_all(lines[i,]$Experiment, "[^[:alnum:]]", "_")
    folder <- paste0(path,"/",name.experiment,"/",name.sample)
    dir.create(folder,showWarnings = F, recursive =T)

    #download files
    for(j in seq_along(files)){
      aux <- paste0(folder,"/",basename(files[j]))
      if(!file.exists(aux)){
        message(paste0("Downloading: ", aux))
        downloader::download(files[j],aux)
      }
    }
  }
  if(length(error)){
    message("=============================")
    message("         WARNING             ")
    message("=============================")
    message("Empty FTP: No files found")
    invisible(apply(as.array(error),1,function(x)  message("Accession: ",x)))
    message("=============================")
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
tcga.download <- function(data=NULL, path=".")
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
