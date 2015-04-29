download <- function(lines){

  encode.url <- "https://www.encodeproject.org/"
  json <- "/?format=JSON"

  for (i in 1:length(lines)){
    if(lines$database[i] == 'encode'){
      id <-  lines$ID[i]
      url <- paste0(encode.url, 'experiments/', id, json)

      print(url)

      item <- rjson::fromJSON (
        RCurl::getURL(url, dirlistonly = TRUE,
                      .opts = list(ssl.verifypeer = FALSE))
      )[['files']]

      for (i in 1:length(item)){
        file <- item[[i]]$href
        link <- paste0(encode.url, file)

        path <- unlist(strsplit(link, "/"))
        path <- path[length(path)]

        print(link)

        downloader::download(link, path, quiet = T)
      }

    }

  }
}
