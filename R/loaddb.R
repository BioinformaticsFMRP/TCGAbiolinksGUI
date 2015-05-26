#' @importFrom rjson fromJSON
load.encode <- function(env) {
    url1 <- c(paste0("https://www.encodeproject.org/search/?type=experiment&",
                     "replicates.library.biosample.donor.organism.scientific_",
                     "name=Homo%20sapiens&limit=all&format=JSON"),
              "Homo sapiens")
    url2 <- c(paste0("https://www.encodeproject.org/search/?type=experiment&",
                     "replicates.library.biosample.donor.organism.scientific_",
                     "name=Mus%20musculus&limit=all&format=JSON"),
              "Mus musculus")
    url3 <- c(paste0("https://www.encodeproject.org/search/?type=experiment&",
                     "replicates.library.biosample.donor.organism.scientific_",
                     "name=Drosophila+melanogaster&limit=all&format=JSON"),
              "Drosophila melanogaster")
    url <- rbind(url1, url2, url3)

    encodeData <- as.data.frame(matrix(ncol = 7))
    colnames(encodeData) <- c("accession", "biosample", "assay",
                              "lab", "target", "description", "organism")

    for (j in 1:dim(url)[1]) {
        message(paste0("Loading ", url[j, 2], "..."))
        data <- fromJSON(getURL(url[j, 1], dirlistonly = TRUE,
                                .opts = list(ssl.verifypeer = FALSE))
        )[["@graph"]]


        pb <- txtProgressBar(style = 3, max = length(data))

        for (i in 1:length(data)) {

            if (is.null(data[[i]]$accession)) {
                accession <- NA
            } else {
                accession <- data[[i]]$accession
            }
            if (is.null(data[[i]]$biosample_term_name)) {
                biosample <- NA
            } else {
                biosample <- data[[i]]$biosample_term_name
            }
            if (is.null(data[[i]]$assay_term_name)) {
                assay <- NA
            } else {
                assay <- data[[i]]$assay_term_name
            }
            if (is.null(data[[i]]$lab$title)) {
                lab <- NA
            } else {
                lab <- data[[i]]$lab$title
            }
            if (is.null(data[[i]]$target$label)) {
                target <- NA
            } else {
                target <- data[[i]]$target$label
            }
            if (is.null(data[[i]]$description)) {
                description <- NA
            } else {
                description <- data[[i]]$description
            }

            encodeData <- rbind(encodeData,
                                c(accession, biosample,
                                  assay, lab, target, description, url[j, 2]))

            setTxtProgressBar(pb, i)
        }
        close(pb)
    }
    encode.db <- encodeData[-1, ]
    assign("encode.db",encode.db, envir = env)
}

load.roadmap <- function(env) {
    dataRoadmap <- read.csv(paste0("http://www.ncbi.nlm.nih.gov/geo/roadmap/",
                                   "epigenomics/?view=samples&sort=",
                                   "acc&mode=csv"),
                            stringsAsFactors = FALSE)
    assign("roadmap.db",dataRoadmap, envir = env)
}

#' @importFrom downloader download
load.biosamples <- function(env) {

    gdocs <- paste0("https://docs.google.com/spreadsheets/d/10GwiiO8A4Ld1h",
                    "4HTGO88oaP7y3sqJHNRiQ_wcnKfXyM/export?format=tsv&id=",
                    "10GwiiO8A4Ld1h4HTGO88oaP7y3sqJHNRiQ_wcnKfXyM&gid=")

    download(paste0(gdocs, "28842439"), destfile = "roadmaptab.tsv")
    download(paste0(gdocs, "514976736"), destfile = "encodetab.tsv")
    download(paste0(gdocs, "0"), destfile = "tcgatab.tsv")

    biosample.encode <- read.delim(file = "encodetab.tsv", sep = "\t")
    biosample.roadmap <- read.delim(file = "roadmaptab.tsv", sep = "\t")
    biosample.tcga <- read.delim(file = "tcgatab.tsv", sep = "\t")

    biosample.encode <- data.frame(lapply(biosample.encode,
                                              as.character),
                                       stringsAsFactors = FALSE)
    biosample.roadmap <- data.frame(lapply(biosample.roadmap,
                                               as.character),
                                        stringsAsFactors = FALSE)
    biosample.tcga <- data.frame(lapply(biosample.tcga, as.character),
                                     stringsAsFactors = FALSE)

    if (file.exists("encodetab.tsv")) {
        file.remove("encodetab.tsv")
    }
    if (file.exists("roadmaptab.tsv")) {
        file.remove("roadmaptab.tsv")
    }
    if (file.exists("tcgatab.tsv")) {
        file.remove("tcgatab.tsv")
    }
    assign("biosample.tcga",biosample.tcga, envir = env)
    assign("biosample.roadmap",biosample.roadmap, envir = env)
    assign("biosample.encode",biosample.encode, envir = env)
}

#' @importFrom downloader download
load.systems <- function(env) {
    gdocs <- paste0("https://docs.google.com/spreadsheets/d/10GwiiO8A4Ld1h",
                    "4HTGO88oaP7y3sqJHNRiQ_wcnKfXyM/export?format=tsv&id=",
                    "10GwiiO8A4Ld1h4HTGO88oaP7y3sqJHNRiQ_wcnKfXyM&gid=")

    download(url = paste0(gdocs, "1660203777"), destfile = "systemstab.tsv")
    systems <- read.delim(file = "systemstab.tsv", sep = "\t")
    systems <- data.frame(lapply(systems, as.character),
                              stringsAsFactors = FALSE)
    if (file.exists("systemstab.tsv")) {
        file.remove("systemstab.tsv")
    }
    assign("systems",systems, envir = env)
}



#' @importFrom downloader download
load.platforms <- function(env) {
    gdocs <- paste0("https://docs.google.com/spreadsheets/d/10GwiiO8A4Ld1h4",
                    "HTGO88oaP7y3sqJHNRiQ_wcnKfXyM/export?format=tsv&id=",
                    "10GwiiO8A4Ld1h4HTGO88oaP7y3sqJHNRiQ_wcnKfXyM&gid=")
    download(paste0(gdocs, "1664832653"), destfile = "platformstab.tsv")
    platforms <- read.delim(file = "platformstab.tsv", sep = "\t")
    platforms <- data.frame(lapply(platforms, as.character),
                                stringsAsFactors = FALSE)
    if (file.exists("platformstab.tsv")) {
        file.remove("platformstab.tsv")
    }
    assign("platforms",platforms, envir = env)
}

# Updates tcga platform and diseases
# @param env package environment
#' @importFrom stringr str_match
#' @importFrom XML readHTMLTable
#' @importFrom downloader download
#' @keywords internal
load.tcga <- function(env) {
    tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"

    # Get platform table
    tcga.query <- "query=Platform"
    next.url <- paste0(tcga.root, tcga.query)
    platform.table <- tcgaGetTable(next.url)
    platform.table <- platform.table[, 1:4]
    platform.table <- platform.table[order(platform.table$name,
                                           decreasing = TRUE),]

    # Get disease table
    tcga.query <- "query=Disease"
    next.url <- paste0(tcga.root, tcga.query)
    disease.table <- tcgaGetTable(next.url)
    disease.table <- disease.table[, 1:3]

    # Get center table
    tcga.query <- "query=Center"
    next.url <- paste0(tcga.root, tcga.query)
    center.table  <- tcgaGetTable(next.url)
    center.table <- center.table[, 1:3]

    if (file.exists("tcga.html")) {
        file.remove("tcga.html")
    }

    assign("platform.table",platform.table, envir = env)
    assign("disease.table",disease.table, envir = env)
    assign("center.table",center.table, envir = env)

    tcga.db <-  createTcgaTable()
    assign("tcga.db",tcga.db, envir = env)
}
