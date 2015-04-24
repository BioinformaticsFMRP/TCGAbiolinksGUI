load.encode <- function(){
  url1 <- c("https://www.encodeproject.org/search/?type=experiment&replicates.library.biosample.donor.organism.scientific_name=Homo%20sapiens&limit=all&format=JSON", "Homo sapiens")
  url2 <- c("https://www.encodeproject.org/search/?type=experiment&replicates.library.biosample.donor.organism.scientific_name=Mus%20musculus&limit=all&format=JSON", "Mus musculus")
  url3 <- c("https://www.encodeproject.org/search/?type=experiment&replicates.library.biosample.donor.organism.scientific_name=Drosophila+melanogaster&limit=all&format=JSON", "Drosophila melanogaster")
  url <- rbind(url1, url2, url3)

  encodeData <- as.data.frame(matrix(ncol = 7))
  colnames(encodeData) <- c("accession", "biosample","assay","lab","target","description", "organism")

  for (j in 1:dim(url)[1]){
    message(paste0('Loading ', url[j,2], "..."))
    data <- rjson::fromJSON (
      RCurl::getURL (url[j,1], dirlistonly = TRUE,
                     .opts = list(ssl.verifypeer = FALSE))
    )[['@graph']]


    pb <- txtProgressBar(style = 3, max = length(data))

    for (i in 1:length(data)){

      if(is.null(data[[i]]$accession)){accession <- NA} else {accession <- data[[i]]$accession}
      if(is.null(data[[i]]$biosample_term_name)){biosample <- NA} else {biosample <- data[[i]]$biosample_term_name}
      if(is.null(data[[i]]$assay_term_name)){assay <- NA} else {assay <- data[[i]]$assay_term_name}
      if(is.null(data[[i]]$lab$title)){lab <- NA} else {lab <- data[[i]]$lab$title}
      if(is.null(data[[i]]$target$label)){target <- NA} else {target <- data[[i]]$target$label}
      if(is.null(data[[i]]$description)){description <- NA} else{description <- data[[i]]$description}

      encodeData <- rbind(encodeData, c(accession, biosample, assay, lab, target, description, url[j,2]))

      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  return (encodeData[-1,])
}

load.roadmap <- function(){
  dataRoadmap <- read.csv("http://www.ncbi.nlm.nih.gov/geo/roadmap/epigenomics/?view=samples&sort=acc&mode=csv",
                          stringsAsFactors = FALSE)
  return (dataRoadmap)
}


load.biosamples <- function(biomicsEnv){
  library(downloader)
  downloader::download(url = "https://docs.google.com/spreadsheets/d/10GwiiO8A4Ld1h4HTGO88oaP7y3sqJHNRiQ_wcnKfXyM/export?format=tsv&id=10GwiiO8A4Ld1h4HTGO88oaP7y3sqJHNRiQ_wcnKfXyM&gid=28842439", destfile = 'roadmaptab.tsv')
  downloader::download(url = "https://docs.google.com/spreadsheets/d/10GwiiO8A4Ld1h4HTGO88oaP7y3sqJHNRiQ_wcnKfXyM/export?format=tsv&id=10GwiiO8A4Ld1h4HTGO88oaP7y3sqJHNRiQ_wcnKfXyM&gid=514976736", destfile = 'encodetab.tsv')
  downloader::download(url = "https://docs.google.com/spreadsheets/d/10GwiiO8A4Ld1h4HTGO88oaP7y3sqJHNRiQ_wcnKfXyM/export?format=tsv&id=10GwiiO8A4Ld1h4HTGO88oaP7y3sqJHNRiQ_wcnKfXyM&gid=0", destfile = 'tcgatab.tsv')

  biomicsEnv$biosample.encode <- read.delim(file = 'encodetab.tsv', sep = '\t', )
  biomicsEnv$biosample.roadmap <- read.delim(file = 'roadmaptab.tsv', sep = '\t')
  biomicsEnv$biosample.tcga <- read.delim(file = 'tcgatab.tsv', sep = '\t')

  biomicsEnv$biosample.encode <- data.frame(lapply(biosample.encode, as.character), stringsAsFactors=FALSE)
  biomicsEnv$biosample.roadmap <- data.frame(lapply(biosample.roadmap, as.character), stringsAsFactors=FALSE)
  biomicsEnv$biosample.tcga <- data.frame(lapply(biosample.tcga, as.character), stringsAsFactors=FALSE)

  if (file.exists('encodetab.tsv')) {file.remove('encodetab.tsv')}
  if (file.exists('roadmaptab.tsv')) {file.remove('roadmaptab.tsv')}
  if (file.exists('tcgatab.tsv')) {file.remove('tcgatab.tsv')}
}

load.systems <- function(biomicsEnv){
  downloader::download(url = "https://docs.google.com/spreadsheets/d/10GwiiO8A4Ld1h4HTGO88oaP7y3sqJHNRiQ_wcnKfXyM/export?format=tsv&id=10GwiiO8A4Ld1h4HTGO88oaP7y3sqJHNRiQ_wcnKfXyM&gid=1660203777", destfile = 'systemstab.tsv')
  biomicsEnv$systems <- read.delim(file = 'systemstab.tsv', sep = '\t', )
  biomicsEnv$systems <- data.frame(lapply(systems, as.character), stringsAsFactors=FALSE)
  if (file.exists('systemstab.tsv')) {file.remove('systemstab.tsv')}
}

load.platforms <- function(biomicsEnv ){
  downloader::download(url = "https://docs.google.com/spreadsheets/d/10GwiiO8A4Ld1h4HTGO88oaP7y3sqJHNRiQ_wcnKfXyM/export?format=tsv&id=10GwiiO8A4Ld1h4HTGO88oaP7y3sqJHNRiQ_wcnKfXyM&gid=1664832653", destfile = 'platformstab.tsv')
  biomicsEnv$platforms <- read.delim(file = 'platformstab.tsv', sep = '\t', )
  biomicsEnv$platforms <- data.frame(lapply(platforms, as.character), stringsAsFactors=FALSE)
  if (file.exists('platformstab.tsv')) {file.remove('platformstab.tsv')}
}
