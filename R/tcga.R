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
#------------------------ Encode search
#' @title TCGA search
#' @description
#'    Searches in the tcga database
#' @param tumor Disease Examples:
#' \tabular{lllll}{
#'OV   \tab BRCA \tab CESC \tab ESCA \tab PCPG\cr
#'LUSC \tab LGG  \tab SKCM \tab KICH \tab CHOL\cr
#'GBM  \tab UCEC \tab PRAD \tab PAAD \tab THYM\cr
#'KIRC \tab THCA \tab SARC \tab LAML \tab TGCT\cr
#'COAD \tab KIRP \tab HNSC \tab ACC  \tab UVM \cr
#'READ \tab BLCA \tab DLBC \tab UCS  \tab FPPP\cr
#'LUAD \tab LIHC \tab STAD \tab MESO \tab CNTL
#'}
#' @param platform Example:
#' \tabular{ll}{
#'CGH-1x1M_G4447A                   \tab IlluminaGA_RNASeqV2                \cr
#'AgilentG4502A_07                  \tab IlluminaGA_mRNA_DGE                \cr
#'Human1MDuo                        \tab HumanMethylation450                \cr
#'HG-CGH-415K_G4124A                \tab IlluminaGA_miRNASeq                \cr
#'HumanHap550                       \tab IlluminaHiSeq_miRNASeq             \cr
#'ABI                               \tab H-miRNA_8x15K                      \cr
#'HG-CGH-244A                       \tab SOLiD_DNASeq                       \cr
#'IlluminaDNAMethylation_OMA003_CPI \tab IlluminaGA_DNASeq_automated        \cr
#'IlluminaDNAMethylation_OMA002_CPI \tab HG-U133_Plus_2                     \cr
#'HuEx-1_0-st-v2                    \tab Mixed_DNASeq                       \cr
#'H-miRNA_8x15Kv2                   \tab IlluminaGA_DNASeq_curated          \cr
#'MDA_RPPA_Core                     \tab IlluminaHiSeq_TotalRNASeqV2        \cr
#'HT_HG-U133A                       \tab IlluminaHiSeq_DNASeq_automated     \cr
#'diagnostic_images                 \tab microsat_i                         \cr
#'IlluminaHiSeq_RNASeq              \tab SOLiD_DNASeq_curated               \cr
#'IlluminaHiSeq_DNASeqC             \tab Mixed_DNASeq_curated               \cr
#'IlluminaGA_RNASeq                 \tab IlluminaGA_DNASeq_Cont_automated   \cr
#'IlluminaGA_DNASeq                 \tab IlluminaHiSeq_WGBS                 \cr
#'pathology_reports                 \tab IlluminaHiSeq_DNASeq_Cont_automated\cr
#'Genome_Wide_SNP_6                 \tab bio                                \cr
#'tissue_images                     \tab Mixed_DNASeq_automated             \cr
#'HumanMethylation27                \tab Mixed_DNASeq_Cont_curated          \cr
#'IlluminaHiSeq_RNASeqV2            \tab Mixed_DNASeq_Cont
#'}
#' @param level "1" "2" "3"
#' @param added.since 04-14-2010
#' @param added.up.to 04-14-2010
#' @param listSample List of samples. Example "TCGA-02-*, c("TCGA-04-06-*","TCGA-04-08-*")
#' @export
#' @import downloader
tcga.search <- function(tumor=NULL,
                      platform=NULL,
                      added.since=NULL,
                      added.up.to=NULL,
                      listSample=NULL,
                      level=NULL
){
  if(!is.null(tumor)){
    if(!(is.element(tolower(tumor),tolower(disease.table$abbreviation)))){
      message("Disease not found. Chosse between:")
      message(paste(disease.table$abbreviation, collapse = " "))
      stop("Invalid tumor")
    }
  }
  if(!is.null(platform)){
    if(!(is.element(tolower(platform),tolower(platform.table$alias)))){
      message("Platform not found. Chosse between:")
      message(paste(platform.table$alias, collapse = " "))
      stop("Invalid platform")
    }
  }

  if(!is.null(level)){
    if(!(is.element(level,c("1","2","3")))){
      message("Levelnot found. Chosse between:'1', '2' or '3'")
      stop("Invalid platform")
    }
  }

  if(!is.null(added.since)){
    d <- try( as.Date( added.since, format= "%m/%d/%Y" ) )
    if( class( d ) == "try-error" || is.na( d ) ) print( "Date format should be mm/dd/YYYY" )
  }
  if(!is.null(added.up.to)){
    d <- try( as.Date( added.up.to, format= "%m/%d/%Y" ) )
    if( class( d ) == "try-error" || is.na( d ) ) print( "Date format should be mm/dd/YYYY" )
  }

  # to be improved
  if(!is.null(listSample)){
    #archives <- c()
    files <- c()
    for(i in seq_along(listSample)){
      # table with barcode id
      # example: query=BiospecimenBarcode[@barcode=TCGA-28-2499*]
      message("Searching for barcode files...")
      message(paste("Barcode:",listSample[i]))
      db <- get.barcode.table(listSample[i])

      # Improvement: using portion analyte in order to select platform
      if(!is.null(platform)){
        analyte <- c()
        pat <- "TCGA-[[:alnum:]]{2}-[[:alnum:]]{4}-[[:alnum:]]{3}-[[:alnum:]]{2}"
        if(is.element(tolower(platform),tolower(dna.plat))){
          analyte = c(analyte,"D")
        }
        if(is.element(tolower(platform),tolower(rna.plat))){
          analyte = c(analyte,"R")
        }
        if(is.element(tolower(platform),tolower(total.rna.plat))){
          analyte = c(analyte,"T")
        }
        if(is.element(tolower(platform),tolower(wgarubcon.plat))){
          analyte = c(analyte,"G")
        }
        if(is.element(tolower(platform),tolower(mirna.plat))){
          analyte = c(analyte,"H")
        }
        if(is.element(tolower(platform),tolower(wgaqiagen2.plat))){
          analyte = c(analyte,"X")
        }
        if(is.element(tolower(platform),tolower(wgaqiagen1.plat))){
          analyte = c(analyte,"W")
        }
        idx <- c()
        for(i in seq_along(analyte)){
          aux <- grep(paste0(pat,"[",analyte[i],"]"),as.character(db$barcode))
          idx <- union(idx,aux)
        }
        print(length(idx))
        db <- db[idx,]
      }
      if(!is.null(platform)){
        pat <- "TCGA-[[:alnum:]]{2}-[[:alnum:]]{4}-[[:alnum:]]{3}-[[:alnum:]]{3}-[[:alnum:]]{4}-"
        idx <- grep(platform,plat.center$platform)
        centers <- as.character(plat.center[idx,"center"])
        idx <- c()
        for(i in seq_along(centers)){
          aux <- grep(paste0(pat,centers[i]),as.character(db$barcode))
          idx <- union(idx,aux)
        }
        print(length(idx))
        db <- db[idx,]
      }

      # get getBcrArchiveCollection table
      for(i in seq_along(db$id)){
        #aux <- getBcrArchiveCollection(db[i,"id"])
        #archives <- rbind(archives,aux)
        aux <- get.samples.files(db[i,"id"])
        aux$barcode <- db[i,"barcode"]
        if(!is.null(aux)){
          files <- rbind(files,aux)
        }
      }
    }
    x <- subset(files, files$isLatest == 1)
    x <- x[!duplicated(x), ]
    x <- x[,order(names(x))]
    if(!is.null(platform)){
      x <- subset(x, grepl(tolower(platform),tolower(x$Platform)))
    }
    if(!is.null(platform)){
      x <- subset(x, tolower(Disease) == tolower(tumor))
    }
    if(!is.null(level)){
      x <- subset(x, grepl(paste0("Level_",level),name))
    }

  }
  else{
    message("CREATING TABLE")
    x <- create.tcga.table (platform=platform,type=level,disease=tumor)
  }
  if(!is.null(added.since)){
    x <- subset(x, as.Date(addedDate) > as.Date(added.since,"%m/%d/%Y"))
  }
  if(!is.null(added.up.to)){
    x <- subset(x, as.Date(addedDate) < as.Date(added.up.to,"%m/%d/%Y"))
  }

  return(x)
}

#@description get the arhive info from tcga api
getArchive <- function(id){
  tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
  tcga.query <- paste0("query=Archive&FileInfo[@id=",id,"]&roleName=archiveCollection")
  url <- paste0(tcga.root,tcga.query)
  db <- tcga.get.table(url)
  return(db)
}

get.samples.files <- function(id){
  archives <- NULL
  tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
  tcga.query <- paste0("query=FileInfo&BiospecimenBarcode[@id=",id,"]&roleName=fileCollection")
  url <- paste0(tcga.root,tcga.query)
  db <- tcga.get.table(url)
  if(!is.null(db)){
    obsolete <- grep("-",db$md5sum)

    if(length(obsolete)>0){db <- db[-obsolete,]}
    files <- unique(db$name)
    #for each file, get the highest ID
    for(i in seq_along(files)){
      same.files <- grep(files[i],db$name)
      latest <- max(as.integer(db[same.files,]$id))
      if(!exists("archives")){
        archives  <- getArchive(latest)
        archives$file <- files[i]

      } else {
        aux  <- getArchive(latest)
        aux$file <- files[i]
        archives <- rbind(archives,aux)
      }
    }
    archives <- archives[,-c(10:15)]
    archives$addedDate <- as.Date(archives$addedDate,"%m-%d-%Y")
    archives <- tcga.db.addCol(archives)
  }
  return(archives)
}

getBcrArchiveCollection <- function (id){
  tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
  tcga.query <- paste0("query=Archive&BiospecimenBarcode[@id=",id,"]&roleName=bcrArchiveCollection")
  url <- paste0(tcga.root,tcga.query)
  db <- tcga.get.table(url)
  db$addedDate <- as.Date(db$addedDate,"%m-%d-%Y")
  db <- tcga.db.addCol(db)
}
get.barcode.table <- function(barcode){
  tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
  tcga.query <- paste0("query=BiospecimenBarcode[@barcode=",barcode,"][@isValid=true]")
  url <- paste0(tcga.root,tcga.query)
  db <- tcga.get.table(url)
}
# This function created tcga.db
# warning to create the db takes almost 3 days
# to update it should take some minutes (function still needed to be done)
#' @import downloader XML plyr stringr
create.tcga.table <- function(disease=NULL, platform=NULL,type=NULL){
  # get all compressed archives
  tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
  regex <- '<table summary="Data Summary".*</a></td></tr></table>'
  message("Downloading TCGA database")
  tcga.query <- paste0("query=Archive[isLatest=1]")

  extra <- ""
  pages <- 32
  if(!is.null(platform)){
    extra <- paste0(extra,"[Platform[@name=",platform,"]]")
    pages <- pages/4
  }
  if(!is.null(disease)){
    extra <- paste0(extra,"[Disease[@abbreviation=",disease,"]]")
    pages <- pages/4
  }
  if(!is.null(type)){
    extra <- paste0(extra,"[ArchiveType[@type=Level_",type,"]]")
    pages <- pages/4
  }
  url <- paste0(tcga.root,tcga.query,extra)
  db <- tcga.get.table(url,pages)

  # remove useless cols
  db <- db[,1:9]

  # remove protected data from the update
  idx <- grep("tcga4yeo",db$deployLocation)
  if(length(idx)>0){
    db <- db[-idx,]
  }
  # transoform date into same format
  db$addedDate <- as.Date(db$addedDate,"%m-%d-%Y")
  db <- tcga.db.addCol(db)
  return(db)
}

# Using the name create two collumns Platform and Disease
tcga.db.addCol <- function(data){
  data$Platform <- ""
  data$Disease <- ""
  diseases <- disease.table$abbreviation
  for(i in seq_along(diseases)){
    idx <- grep(diseases[i],data$baseName)
    if(length(idx >0)){
      data[idx,]$Disease <- diseases[i]
    }
  }

  for(i in seq_along(platform.table$name)){
    idx <- grep(platform.table[i,]$name,data$baseName)
    if(length(idx)>0){
      data[idx,]$Platform <- platform.table[i,]$name
    }
  }
  return (data)
}

# A function to get a tcga table from api
# input: url (max for progress bar)
# return: table
tcga.get.table <- function(url,max=0){
  db <- NULL
  next.url <- url
  regex <- '<table summary="Data Summary".*</a></td></tr></table>'
  next.regex <-"http://tcga-data.nci.nih.gov/tcgadccws/GetHTML.*Next"

  if(max>0){
    pb <- txtProgressBar(min = 0, max = max, style = 3)
    i <- 0
  }
  #print(url)
  while(!is.na(next.url)){
    downloader::download(next.url,"tcga.html",quiet=T)
    html <- readLines("tcga.html")
    match <- str_match(html,regex)
    idx <- which(!is.na(match))
    if(length(idx)==0){
      next.url <- NA
      break
    }
    table <- readHTMLTable(toString(match[idx,]),
                           header = T,
                           stringsAsFactors = FALSE)$'NULL'
    colnames(table) <-table[1,]
    table <- table[-1,]

    if(exists("db")) {
      db <- rbind(db,table)
    } else {
      db <- table
    }
    # get next table
    next.url <- str_match(html,next.regex)[6,]
    next.url <- gsub("amp;","",gsub('\">Next',"",next.url))

    if(max>0){
      # update progress bar
      i <- i + 1
      setTxtProgressBar(pb, i)
    }
  }
  if(max>0){
    setTxtProgressBar(pb, max)
    close(pb)
  }
  return (db)
}

