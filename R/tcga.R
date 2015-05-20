#' @importFrom gdata trim
getBarCode <- function(TCGAList, filter) {
    list <- list[, 2]
    files <- NULL
    for (i in unlist(filter)) {
        for (j in list) {
            if (grepl(i, j))
                files <- paste0(files, " ", j)
        }
    }
    files <- gsub(" ", ",", trim(files))
}
# ------------------------ Encode search
#' @title TCGAQuery
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
#'CGH- 1x1M_G4447A                   \tab IlluminaGA_RNASeqV2   \cr
#'AgilentG4502A_07                  \tab IlluminaGA_mRNA_DGE    \cr
#'Human1MDuo                        \tab HumanMethylation450    \cr
#'HG-CGH-415K_G4124A                \tab IlluminaGA_miRNASeq    \cr
#'HumanHap550                       \tab IlluminaHiSeq_miRNASeq \cr
#'ABI                               \tab H-miRNA_8x15K  \cr
#'HG-CGH-244A                       \tab SOLiD_DNASeq                       \cr
#'IlluminaDNAMethylation_OMA003_CPI \tab IlluminaGA_DNASeq_automated        \cr
#'IlluminaDNAMethylation_OMA002_CPI \tab HG-U133_Plus_2                 \cr
#'HuEx- 1_0-st-v2 \tab Mixed_DNASeq                  \cr
#'H-miRNA_8x15Kv2 \tab IlluminaGA_DNASeq_curated      \cr
#'MDA_RPPA_Core   \tab IlluminaHiSeq_TotalRNASeqV2    \cr
#'HT_HG-U133A     \tab IlluminaHiSeq_DNASeq_automated \cr
#'diagnostic_images                 \tab microsat_i                     \cr
#'IlluminaHiSeq_RNASeq              \tab SOLiD_DNASeq_curated           \cr
#'IlluminaHiSeq_DNASeqC             \tab Mixed_DNASeq_curated           \cr
#'IlluminaGA_RNASeq                 \tab IlluminaGA_DNASeq_Cont_automated   \cr
#'IlluminaGA_DNASeq                 \tab IlluminaHiSeq_WGBS             \cr
#'pathology_reports                 \tab IlluminaHiSeq_DNASeq_Cont_automated\cr
#'Genome_Wide_SNP_6                 \tab bio                            \cr
#'tissue_images                     \tab Mixed_DNASeq_automated         \cr
#'HumanMethylation27                \tab Mixed_DNASeq_Cont_curated      \cr
#'IlluminaHiSeq_RNASeqV2            \tab Mixed_DNASeq_Cont
#'}
#' @param level '1' '2' '3'
#' @param added.since 04- 14-2010
#' @param added.up.to 04- 14-2010
#' @param center center name
#' @param samples List of samples. Ex:c('TCGA-04-06-*','TCGA-04-08-*')
#' @example inst/examples/tcgaSearch.R
#' @export
#' @importFrom downloader download
#' @importFrom knitr kable
#' @return A dataframe with the results of the query
#'        (lastest version of the files)
tcgaSearch <- function(tumor = NULL, platform = NULL, added.since = NULL,
                      added.up.to = NULL, samples = NULL, center = NULL,
                      level = NULL) {
    disease.table   <- get("disease.table")
    platform.table  <- get("platform.table")
    center.table  <- get("center.table")
    db <-  get("tcga.db")

    if (!is.null(tumor)) {
        for (j in seq_along(tumor)) {
            if (!(is.element(tolower(tumor[j]),
                             tolower(disease.table$abbreviation)))) {
                suppressWarnings(
                    df <- as.data.frame(matrix(sort(unique(disease.table$abbreviation)),
                                               ncol = 8))
                )
                print(kable(df, col.names = NULL, format = "pandoc",
                            caption = "TCGA tumors"))
                cat("=======================================================\n")
                cat("ERROR: Disease not found. Select from the table above.\n")
                cat("=======================================================\n")
                return(NULL)
            }
        }
    }

    if (!is.null(platform)) {
        for (j in seq_along(platform)) {
            if (!(is.element(tolower(platform[j]), tolower(platform.table$name)))) {
                suppressWarnings(
                    df <- as.data.frame(matrix(sort(unique(platform.table$name)),
                                               ncol = 3))
                )
                print(kable(df, col.names = NULL, format = "pandoc",
                            caption = "TCGA Platforms"))
                cat("=======================================================\n")
                cat("ERROR: Platform not found. Select from the table above.\n")
                cat("=======================================================\n")
                return(NULL)
            }
        }
    }

    if (!is.null(center)) {
        if (!(is.element(tolower(center), tolower(center.table$name)))) {
            suppressWarnings(
                df <- as.data.frame(matrix(sort(unique(center.table$name)),
                                           ncol = 3))
            )
            print(kable(df, col.names = NULL, format = "pandoc",
                        caption = "TCGA Centers"))
            cat("=======================================================\n")
            cat("ERROR: Center not found. Select from the table above.\n")
            cat("=======================================================\n")
            return(NULL)
        }
    }
    if (!is.null(level)) {
        if (!(is.element(level, c("1", "2", "3","mage-tab")))) {
            message("Level not found. Chosse between:'1', '2','3','mage-tab'")
            stop("Invalid platform")
        }
    }

    if (!is.null(added.since)) {
        d <- try(as.Date(added.since, format = "%Y/%m/%d"))
        if (class(d) == "try-error" || is.na(d)) {
            print("Date format should be YYYY-mm-dd")
        }
    }
    if (!is.null(added.up.to)) {
        d <- try(as.Date(added.up.to, format = "%Y/%m/%d"))
        if (class(d) == "try-error" || is.na(d)) {
            print("Date format should be YYYY-mm-dd")
        }
    }

    if(!is.null(tumor)){
        id <- sapply(tumor, function(x){
            grepl(x, db$Disease, ignore.case = TRUE)
        })
        id <- apply(id, 1,any)
        db <-  db[id,]
    }
    if(!is.null(platform)){
        id <- sapply(platform, function(x){
            grepl(x, db$Platform, ignore.case = TRUE)
        })
        id <- apply(id, 1,any)
        db <-  db[id,]
    }
    if(!is.null(center)){
        id <- sapply(center, function(x){
            grepl(x, db$Center, ignore.case = TRUE)
        })
        id <- apply(id, 1,any)
        db <-  db[id,]
    }
    if(!is.null(level)){
        id <- grep(paste0("Level_", level), db$name)
        if(length(id) > 0){
            db <-  db[id,]
        }
    }

    if (!is.null(added.since)) {
        db <- subset(db, as.Date(db$addedDate) > as.Date(added.since,
                                                         "%m/%d/%Y"))
    }
    if (!is.null(added.up.to)) {
        db <- subset(db, as.Date(db$addedDate) < as.Date(added.up.to,
                                                         "%m/%d/%Y"))
    }

    # to be improved
    idx <- c()
    if(!is.null(samples)){
        for(i in seq_along(samples)){
            aux <- grep(samples[i],db$barcode)
            idx <- union(idx, aux)
        }
        db <- db[idx,]
    }
    return(db)
}
# @description get the arhive info from tcga api
getArchive <- function(id) {
    root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
    query <- paste0("query=Archive&FileInfo[@id=", id,
                    "]&roleName=archiveCollection")
    url <- paste0(root, query)
    db <- tcgaGetTable(url)
    return(db)
}

getSamplesFiles <- function(id) {
    archives <- NULL
    root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
    query <- paste0("query=FileInfo&BiospecimenBarcode[@id=",
                    id, "]&roleName=fileCollection")
    url <- paste0(root, query)
    db <- tcgaGetTable(url)
    if (!is.null(db)) {
        obsolete <- grep("-", db$md5sum)

        if (length(obsolete) > 0) {
            db <- db[-obsolete, ]
        }
        files <- unique(db$name)
        # for each file, get the highest ID
        for (i in seq_along(files)) {
            same.files <- grep(files[i], db$name)
            latest <- max(as.integer(db[same.files, ]$id))
            if (!exists("archives")) {
                archives <- getArchive(latest)
                archives$file <- files[i]

            } else {
                aux <- getArchive(latest)
                aux$file <- files[i]
                archives <- rbind(archives, aux)
            }
        }
        archives <- archives[, -c(10:15)]
        archives$addedDate <- as.Date(archives$addedDate, "%m-%d-%Y")
        archives <- tcgaDbAddCol(archives)
    }
    return(archives)
}

getBcrArchiveCollection <- function(id) {
    tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
    tcga.query <- paste0("query=Archive&BiospecimenBarcode[@id=",
                        id, "]&roleName=bcrArchiveCollection")
    url <- paste0(tcga.root, tcga.query)
    db <- tcgaGetTable(url)
    db$addedDate <- as.Date(db$addedDate, "%m-%d-%Y")
    db <- tcgaDbAddCol(db)
}
getBarcodeTable <- function(barcode) {
    tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
    tcga.query <- paste0("query=BiospecimenBarcode[@barcode=",
                        barcode, "][@isValid=true]")
    url <- paste0(tcga.root, tcga.query)
    db <- tcgaGetTable(url)
    return (db)
}
# This function created tcga.db warning to create the db
# takes almost 3 days to update it should take some minutes
# (function still needed to be done)
#' @importFrom downloader download
createTcgaTable <- function(disease = NULL, platform = NULL, type = NULL) {
    # get all compressed archives
    tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
    message("Downloading TCGA database")
    tcga.query <- paste0("query=Archive[isLatest=1]")

    extra <- ""
    pages <- 32
    if (!is.null(platform)) {
        extra <- paste0(extra, "[Platform[@name=", platform,
                        "]]")
        pages <- pages / 4
    }
    if (!is.null(disease)) {
        extra <- paste0(extra, "[Disease[@abbreviation=", disease,
                        "]]")
        pages <- pages / 4
    }
    if (!is.null(type)) {
        extra <- paste0(extra, "[ArchiveType[@type=Level_", type,
                        "]]")
        pages <- pages / 4
    }
    url <- paste0(tcga.root, tcga.query, extra)
    db <- tcgaGetTable(url, pages)

    # remove useless cols
    db <- db[, 1:9]

    # remove protected data from the update
    idx <- grep("tcga4yeo", db$deployLocation)
    if (length(idx) > 0) {
        db <- db[-idx, ]
    }
    # transoform date into same format
    db$addedDate <- as.Date(db$addedDate, "%m-%d-%Y")
    db <- tcgaDbAddCol(db)
    return(db)
}

# Using the name create two collumns Platform and Disease
tcgaDbAddCol <- function(data) {
    disease.table  <- get("disease.table")
    platform.table  <- get("platform.table")
    data$Platform <- ""
    data$Disease <- ""
    diseases <- disease.table$abbreviation
    for (i in seq_along(diseases)) {
        idx <- grep(diseases[i], data$baseName)
        if (length(idx > 0)) {
            data[idx, ]$Disease <- diseases[i]
        }
    }

    for (i in seq_along(platform.table$name)) {
        idx <- grep(platform.table[i, ]$name, data$baseName)
        if (length(idx) > 0) {
            data[idx, ]$Platform <- platform.table[i, ]$name
        }
    }
    return(data)
}

# A function to get a tcga table from api input: url (max for
# progress bar) return: table
#' @importFrom stringr str_match
#' @importFrom downloader download
tcgaGetTable <- function(url, max = 0) {
    db <- NULL
    next.url <- url
    regex <- "<table summary=\"Data Summary\".*</a></td></tr></table>"
    next.regex <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML.*Next"

    if (max > 0) {
        pb <- txtProgressBar(min = 0, max = max, style = 3)
        i <- 0
    }
    # print(url)
    while (!is.na(next.url)) {
        download(next.url, "tcga.html", quiet = TRUE)
        html <- readLines("tcga.html")
        match <- str_match(html, regex)
        idx <- which(!is.na(match))
        if (length(idx) == 0) {
            next.url <- NA
            break
        }
        table <- readHTMLTable(toString(match[idx, ]), header = TRUE,
                               stringsAsFactors = FALSE)$"NULL"
        colnames(table) <- table[1, ]
        table <- table[-1, ]

        if (exists("db")) {
            db <- rbind(db, table)
        } else {
            db <- table
        }
        # get next table
        next.url <- str_match(html, next.regex)[6, ]
        next.url <- gsub("amp;", "", gsub("\">Next", "", next.url))

        if (max > 0) {
            # update progress bar
            i <- i + 1
            setTxtProgressBar(pb, i)
        }
    }
    if (max > 0) {
        setTxtProgressBar(pb, max)
        close(pb)
    }
    return(db)
}
