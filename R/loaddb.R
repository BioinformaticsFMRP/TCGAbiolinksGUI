#' @importFrom rjson fromJSON
load.encode <- function() {
    url1 <- c(paste0("https://www.encodeproject.org/search/?type=experiment&",
                     "replicates.library.biosample.donor.organism.scientific_",
                     "name=Homo%20sapiens&limit=all&format=JSON"), "Homo sapiens")
    url2 <- c(paste0("https://www.encodeproject.org/search/?type=experiment&",
                     "replicates.library.biosample.donor.organism.scientific_",
                     "name=Mus%20musculus&limit=all&format=JSON"), "Mus musculus")
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
        data <- rjson::fromJSON(RCurl::getURL(url[j, 1], dirlistonly = TRUE,
                                              .opts = list(ssl.verifypeer = FALSE)))[["@graph"]]


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

            encodeData <- rbind(encodeData, c(accession, biosample,
                                              assay, lab, target, description, url[j, 2]))

            setTxtProgressBar(pb, i)
        }
        close(pb)
    }
    return(encodeData[-1, ])
}

load.roadmap <- function() {
    dataRoadmap <- read.csv(paste0("http://www.ncbi.nlm.nih.gov/geo/roadmap/",
                                   "epigenomics/?view=samples&sort=acc&mode=csv"), stringsAsFactors = FALSE)
    return(dataRoadmap)
}

#' @importFrom downloader download
load.biosamples <- function(env) {
    gdocs <- paste0("https://docs.google.com/spreadsheets/d/10GwiiO8A4Ld1h4HTGO",
                    "88oaP7y3sqJHNRiQ_wcnKfXyM/export?format=tsv&id=10GwiiO8A4L",
                    "d1h4HTGO88oaP7y3sqJHNRiQ_wcnKfXyM&gid=")
    download(paste0(gdocs, "28842439"), destfile = "roadmaptab.tsv")
    download(paste0(gdocs, "514976736"), destfile = "encodetab.tsv")
    download(paste0(gdocs, "0"), destfile = "tcgatab.tsv")

    env$biosample.encode <- read.delim(file = "encodetab.tsv",
                                       sep = "\t", )
    env$biosample.roadmap <- read.delim(file = "roadmaptab.tsv",
                                        sep = "\t")
    env$biosample.tcga <- read.delim(file = "tcgatab.tsv", sep = "\t")

    env$biosample.encode <- data.frame(lapply(biosample.encode,
                                              as.character), stringsAsFactors = FALSE)
    env$biosample.roadmap <- data.frame(lapply(biosample.roadmap,
                                               as.character), stringsAsFactors = FALSE)
    env$biosample.tcga <- data.frame(lapply(biosample.tcga, as.character),
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
}

#' @importFrom downloader download
load.systems <- function(env) {
    gdocs <- paste0("https://docs.google.com/spreadsheets/d/10GwiiO8A4Ld1h4HTGO",
                    "88oaP7y3sqJHNRiQ_wcnKfXyM/export?format=tsv&id=10GwiiO8A4L",
                    "d1h4HTGO88oaP7y3sqJHNRiQ_wcnKfXyM&gid=")

    download(url = paste0(gdocs, "1660203777"), destfile = "systemstab.tsv")
    env$systems <- read.delim(file = "systemstab.tsv", sep = "\t")
    env$systems <- data.frame(lapply(systems, as.character),
                              stringsAsFactors = FALSE)
    if (file.exists("systemstab.tsv")) {
        file.remove("systemstab.tsv")
    }
}

#' @importFrom downloader download
load.platforms <- function(env) {
    gdocs <- paste0("https://docs.google.com/spreadsheets/d/10GwiiO8A4Ld1h4HTGO",
                    "88oaP7y3sqJHNRiQ_wcnKfXyM/export?format=tsv&id=10GwiiO8A4L",
                    "d1h4HTGO88oaP7y3sqJHNRiQ_wcnKfXyM&gid=")
    download(paste0(gdocs, "1664832653"), destfile = "platformstab.tsv")
    env$platforms <- read.delim(file = "platformstab.tsv", sep = "\t")
    env$platforms <- data.frame(lapply(platforms, as.character),
                                stringsAsFactors = FALSE)
    if (file.exists("platformstab.tsv")) {
        file.remove("platformstab.tsv")
    }
}

#' Updates tcga platform and diseases
#' @param env package environment
#' @importFrom stringr str_match
#' @importFrom XML readHTMLTable
#' @importFrom downloader download
load.tcga <- function(env) {
    tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
    tcga.query <- "query=Platform"
    next.url <- paste0(tcga.root, tcga.query)
    download(next.url, "tcga.html", quiet = TRUE)
    regex <- "<table summary=\"Data Summary\".*</a></td></tr></table>"
    html <- readLines("tcga.html")
    platform.table <- readHTMLTable(toString(str_match(html,
                                                       regex)[6, ]), header = TRUE, stringsAsFactors = FALSE)$"NULL"
    colnames(platform.table) <- platform.table[1, ]
    env$platform.table <- platform.table[-1, 1:4]

    tcga.query <- "query=Disease"
    next.url <- paste0(tcga.root, tcga.query)
    download(next.url, "tcga.html", quiet = TRUE)
    regex <- "<table summary=\"Data Summary\".*</a></td></tr></table>"
    html <- readLines("tcga.html")
    match <- str_match(html, regex)
    idx <- which(!is.na(match))
    if (length(idx) > 0) {
        disease.table <- readHTMLTable(toString(match[idx, ]),
                                       header = TRUE, stringsAsFactors = FALSE)$"NULL"
        colnames(disease.table) <- disease.table[1, ]
        env$disease.table <- disease.table[-1, 1:4]
    }
    if (file.exists("tcga.html")) {
        file.remove("tcga.html")
    }

    env$tcga.db <- tcga.search()

}
