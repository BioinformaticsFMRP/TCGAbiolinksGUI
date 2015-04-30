#' @export
biOMICs.download <- function(lines){

  encode.url <- "https://www.encodeproject.org/"
  tcga.url <- "https://tcga-data.nci.nih.gov"
  json <- "/?format=JSON"

  for (i in 1:dim(lines)[1]){
    #download Encode
    if(lines$database[i] == 'encode'){
      id <-  lines$ID[i]
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

    #download Roadmap
    if(lines$database[i] == 'roadmap'){
      id <-  lines$ID[i]
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

    #download TCGA
    if(lines$database[i] == 'tcga'){
      id <-  lines$ID[i]
      path <- tcga.db[tcga.db$name == id,]$deployLocation

      link <- paste0(tcga.url, path)

      fileout <- unlist(strsplit(link, "/"))
      fileout <- fileout[length(fileout)]

      downloader::download(link, fileout)
    }
  }
}
