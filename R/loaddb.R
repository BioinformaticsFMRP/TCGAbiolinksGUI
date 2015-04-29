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

#' @import XML  stringr
load.tcga <- function(biomicsEnv){
  downloader::download(url = "https://docs.google.com/spreadsheets/d/10GwiiO8A4Ld1h4HTGO88oaP7y3sqJHNRiQ_wcnKfXyM/export?format=tsv&id=10GwiiO8A4Ld1h4HTGO88oaP7y3sqJHNRiQ_wcnKfXyM&gid=1340244739", destfile = 'tcga.tsv')
  tcga.db <- read.delim(file = 'tcga.tsv', sep = '\t' )
  biomicsEnv$tcga.db <- data.frame(lapply(tcga.db, as.character), stringsAsFactors=FALSE)
  if (file.exists('tcga.tsv')) {file.remove('tcga.tsv')}

  tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
  tcga.query <- "query=Platform"
  next.url <- paste0(tcga.root,tcga.query)
  downloader::download(next.url,"tcga.html",quiet =T)
  regex <- '<table summary="Data Summary".*</a></td></tr></table>'
  html <- readLines("tcga.html")
  platform.table <- readHTMLTable(toString(str_match(html,regex)[6,]),
                                  header = T,
                                  stringsAsFactors = FALSE)$'NULL'
  colnames(platform.table) <- platform.table[1,]
  biomicsEnv$platform.table <- platform.table[-1,1:4]
  if (file.exists('tcga.html')) {file.remove('tcga.html')}

}

#' @import downloader XML plyr stringr
create.tcga.table <- function(biomicsEnv){
  # get all compressed archives
  tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
  regex <- '<table summary="Data Summary".*</a></td></tr></table>'
  message("Downloading TCGA database")
  tcga.query <- paste0("query=Archive[isLatest=1]")
  next.url <- paste0(tcga.root,tcga.query)
  pb <- txtProgressBar(min = 0, max = 32, style = 3)
  i <- 0
  while(!is.na(next.url)){
    downloader::download(next.url,"tcga.html",quiet=T)
    html <- readLines("tcga.html")
    table <- readHTMLTable(toString(str_match(html,regex)[6,]),
                           header = T,
                           stringsAsFactors = FALSE)$'NULL'
    colnames(table) <-table[1,]
    table <- table[-1,1:9]

    if(exists("tcga.db")) {
      tcga.db <- rbind(tcga.db,table)
    } else {
      tcga.db <- table
    }
    # get next table
    next.regex <-"http://tcga-data.nci.nih.gov/tcgadccws/GetHTML.*Next"
    next.url <- str_match(html,next.regex)[6,]
    next.url <- gsub("amp;","",gsub('\">Next',"",next.url))

    # update progress bar
    i <- i + 1
    setTxtProgressBar(pb, i)
  }
  close(pb)
  biomicsEnv$tcga.db <- tcga.db
}


# save as: list (list(barcodes)) or collapse (paste0(barcode, collapse = ","))?
load.tcga.barcode <- function(){
  # todo considere next page =/
  # comparar com sdrf
  # for each tcga.db id get barcodes
  message("Downloading TCGA barcodes")
  start.time <- Sys.time()
  pb <- txtProgressBar(min = 0, max = nrow(tcga.db), style = 3)
  for(j in 1:nrow(tcga.db)){
  #for(j in c(5)){
    tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?query="
    tcga.query <- paste0("FileInfo&Archive[@id=",tcga.db[j,'id'],"]&roleName=fileCollection")
    # metadata files ignore barcodes
    if(length(grep(".*(aux|mage-tab).*",tcga.db[j,'name']))>0){
      print(paste("Continued for:",j))
      next
    }
    # search if a file with same serial index and revision have barcodes
    # if it has we should copy it!
    rev <- tcga.db[j,"revision"]
    index <- tcga.db[j,"serialIndex"]
    bname<- tcga.db[j,"baseName"]
    aux <- subset(tcga.db,
                  revision == rev & serialIndex == index & baseName==bname & deployStatus!="Available",
                  select = c(deployLocation,deployStatus)
                  )

    if(nrow(aux)>0){
      print("Found barcode before")
      tcga.db[j,"deployStatus"] <- aux[1,"deployStatus"]
      next
    }
    next.url <- paste0(tcga.root,tcga.query)
    while(!is.na(next.url)){
      print(next.url)
      downloader::download(next.url,"tcga.html",quiet =T)
      regex <- '<table summary="Data Summary".*</a></td></tr></table>'
      html <- readLines("tcga.html")
      files <- readHTMLTable(toString(str_match(html,regex)[6,]),
                             header = T,
                             stringsAsFactors = FALSE)$'NULL'
      colnames(files) <- files[1,]
      files <- files[-1,1:5]
      idx <- grep("(README|CHANGES|DESCRIPTION|MANIFEST).*",files$name)
      if(length(idx > 0)){files <- files[-idx,]}
      if(nrow(files) == 0) {
        next.url <- NA
        if(!exists("all.barcode")){
          all.barcode <- ""
        }
        next
      }
      pat <- "*(TCGA)-([A-Z0-9]{2})-([A-Z0-9]{4})-(0[1-9]|[1-2][0-9])([A-Z])-(0[1-9]|[1-9][0-9])([DGHRTWX])-([A-Z0-9]{4})-([A-Z0-9]{2})*"
      barcode <- str_match(files$name,pat)[,1]
      #message("Found in name")
      #print(barcode)
      if(all(is.na(barcode))){
        for(i in 1:nrow( files )){
          tcga.query <- paste0("BiospecimenBarcode&FileInfo[@id=",files[i,"id"],"]&roleName=biospecimenBarcodeCollection")
          next.url <- paste0(tcga.root,tcga.query)
          downloader::download(next.url,"tcga.html",quiet =T)
          regex <- '<table summary="Data Summary".*</a></td></tr></table>'
          barcode.html <- readLines("tcga.html")
          table <- str_match(barcode.html,regex)
          idx <- which(!is.na(table))
          if(length(idx>0)){
            barcode.table <- readHTMLTable(toString(table[idx,]),
                                           header = T,
                                           stringsAsFactors = FALSE)$'NULL'
            colnames(barcode.table) <- barcode.table[1,]
            barcode.table <- barcode.table[-1,1:8]
            #print(barcode.table$barcode)
            files[[i,5]] <- list(barcode.table$barcode)
          }
        }
        barcode <- unique(unlist(files[,5]))
      }
      # get next table
      next.regex <-"http://tcga-data.nci.nih.gov/tcgadccws/GetHTML.*Next"
      next.url <- str_match(html,next.regex)[6,]
      next.url <- gsub("amp;","",gsub('\">Next',"",next.url))
      if(exists("all.barcode")){
        all.barcode <- c(all.barcode,barcode)
      } else {
        all.barcode <- barcode
      }
    }
    #tcga.db[[j,"deployStatus"]] <- list(barcode)
    tcga.db[j,"deployStatus"] <- paste0(all.barcode, collapse = ",")
    rm(all.barcode)
    setTxtProgressBar(pb, j)
    print(Sys.time()-start.time)
  }
  close(pb)
}


# Using the name create two collumns Platform and Disease
tcga.db.addCol <- function(x){
  tcga.db$Platform <- ""
  tcga.db$Disease <- ""
  diseases <- sapply(strsplit(biosample.tcga$biosample,split = " - "),function(x){x[1]})
  for(i in seq_along(diseases)){
    idx <- grep(diseases[i],tcga.db$baseName)
    tcga.db[idx,]$Disease <- diseases[i]
  }

  for(i in seq_along(platform.table$name)){
    idx <- grep(platform.table[i,]$name,tcga.db$baseName)
    if(length(idx)>0){
      tcga.db[idx,]$Platform <- platform.table[i,]$name
    }
  }
}

