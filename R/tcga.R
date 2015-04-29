getBarCode <- function(TCGAList, filterList){
  TCGAList <- TCGAList[,2]
  files <- NULL
  for (i in unlist (filterList)) {
    for (j in TCGAList) {
      if (grepl(i,j)) files <- paste0 (files," ", j)
    }
  }
  files <- gsub(" ", ",", gdata::trim(files))
}


tcga.query <- function(tumor="all",
                       platform="all",
                       added.since=NULL,
                       added.up.to=NULL,
                       listSample=NULL,
                       level = "all"
){
  x <- tcga.db
  if(tumor!="all"){
    x <- subset(x, tolower(Disease) == tolower(tumor))
  }
  if(platform!="all"){
    x <- subset(x, tolower(Platform) == tolower(platform))
  }
  if(level!="all"){
    x <- subset(x, grepl(paste0("Level_",level),name) )
  }
  if(!is.null(added.since)){
    x <- subset(x, as.Date(addedDate,"%m/%d/%Y") > as.Date(added.since,"%m/%d/%Y"))
  }
  if(!is.null(added.up.to)){
    x <- subset(x, as.Date(addedDate,"%m/%d/%Y") < as.Date(added.up.to,"%m/%d/%Y"))
  }

  # to be improved
  if(!is.null(listSample)){
    for(i in seq_along(listSample)){
      aux <- grep(listSample[i],tcga.db$deployStatus)
      if(!exists("idx")){
        idx <- aux
      } else {
        idx <- union(idx, aux)
      }
    }
    x <- x[idx,]
  }
  return(x)
}

tcga.download<- function(data=NULL,path=".")
{
  dir.create(path,showWarnings = F)
  root <- "https://tcga-data.nci.nih.gov"
  for(i in 1:nrow(data)){
    file <- paste0(path,"/",basename(data[i,"deployLocation"]))
    message(paste0("Downloading:", basename(data[i,"deployLocation"])))
    downloader::download(paste0(root,x[i,"deployLocation"]),file)
    untar(file, exdir = path)
  }
}
