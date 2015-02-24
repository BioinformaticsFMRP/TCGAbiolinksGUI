.getBarCode <- function(TCGAList, filterList){
  TCGAList <- TCGAList[,2]
  files <- NULL
  for (i in unlist (filterList)) {
    for (j in TCGAList) {
      if (grepl(i,j)) files <- paste0 (files," ", j)
    }
  }
  files <- gsub(" ", ",", trim(files))
}
