#' @title biOmicsDownload
#' @description Download data selected using the biOmics.search function
#' @param lines biOmicsSearch output
#' @param enc.file.type Extension to be downloaed from encode database
#' @param rmap.file.type Extension to be downloaed from roadmap database
#' @seealso biOmicsSearch
#' @export
biOmicsDownload <- function(lines, enc.file.type = NULL,
                            rmap.file.type = NULL) {
    encode.db  <- getOption("encode.db")
    tcga.db    <- getOption("tcga.db")
    roadmap.db <- getOption("roadmap.db")

    with(lines,{
        encode.lines <- subset(lines, lines$database == "encode")
        rmap.lines <- subset(lines, lines$database == "roadmap")
        tcga.lines <- subset(lines, lines$database == "tcga")
    })
    # -------------- download Encode
    if (dim(encode.lines)[1] > 0) {
        encode.lines <- subset(encode.db,
                                encode.db$accession == encode.lines$ID)
        encodeDownload(encode.lines, enc.file.type)
    }

    # -------------- download ROADMAP
    if (dim(rmap.lines)[1] > 0) {
        rmap.lines <- subset(roadmap.db,
                            roadmap.db$X..GEO.Accession == rmap.lines$ID)
        roadmapDownload(rmap.lines, rmap.file.type)
    }
    # ---------------- download TCGA TODO: add filters, folder to
    # save
    if (dim(tcga.lines)[1] > 0) {
        tcga.lines <- tcga.db[tcga.db$name == tcga.lines$ID,]
        tcgaDownload(tcga.lines, path = "TCGA")
    }
}

#' @title Download encode data
#' @description Download encode data selected using the encodeSearch
#' @param lines encode.search output
#' @param path Folder to save the file
#' @param type extesion of files to be downloaded
#' @export
#' @importFrom downloader download
#' @importFrom RCurl getURL
#' @importFrom rjson fromJSON
encodeDownload <- function(lines, type = NULL, path = ".") {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)

    encode.url <- "https://www.encodeproject.org/"
    json <- "/?format=JSON"

    for (i in 1:dim(lines)[1]) {
        id <- lines$accession[i]
        url <- paste0(encode.url, "experiments/", id, json)

        # get list of files
        item <- fromJSON(getURL(url, dirlistonly = TRUE,
                                .opts = list(ssl.verifypeer = FALSE)
                                )
                         )[["files"]]

        files <- sapply(item, function(x) {
            x$href
        })

        if (!is.null(type)) {
            idx <- unlist(lapply(type, function(x) {
                grep(x, files)
            }))
            files <- files[idx]
        }
        # download files
        for (j in seq_along(files)) {
            link <- paste0(encode.url, files[j])
            fileout <- file.path(path, id, basename(link))
            download(link, fileout)
        }
    }

}

#' @title Download roadmap data
#' @description Download roadmap data selected using the roadmap.search
#' @param lines roadmap.search output
#' @param path Folder to save the file
#' @param type extesion of files to be downloaded
#' @export
#' @importFrom downloader download
#' @importFrom stringr str_replace_all
roadmapDownload <- function(lines, type = NULL, path = ".") {

    error <- c()
    for (i in 1:dim(lines)[1]) {
        sample <- str_replace_all(lines[i, ]$Sample.Name, "[^[:alnum:]]","_")
        expr <- str_replace_all(lines[i, ]$Experiment, "[^[:alnum:]]","_")

        folder <- paste0(path, "/", expr, "/", sample)
        status <- dir.create(folder, showWarnings = FALSE, recursive = TRUE)
        if (!status) {
            message("Downloaded", i)
            next
        }

        url <- lines[i, ]$GEO.FTP
        if (nchar(url) == 0) {
            error <- c(error, lines[i, ]$X..GEO.Accession)
            next
        }
        filenames <- getFileNames(url)
        if (!is.null(type)) {
            idx <- unlist(lapply(type, function(x) {
                grep(x, filenames)
            }))
            filenames <- filenames[idx]
        }
        files <- paste0(url, filenames)

        # download files
        for (j in seq_along(files)) {
            aux <- paste0(folder, "/", basename(files[j]))
            if (!file.exists(aux)) {
                message(paste0("Downloading: ", aux))
                download(files[j], aux)
            }
        }
    }
    if (length(error)) {
        message("=============================")
        message("         WARNING             ")
        message("=============================")
        message("Empty FTP: No files found")
        invisible(apply(as.array(error), 1,
                        function(x) message("Accession: ",x))
                  )
        message("=============================")
    }
}
#' @title TCGA Download
#' @description Download data previously selected using the TCGASeach
#' @param data TCGASearch output
#' @param path location of the final data saving
#' @seealso TCGASearch
#' @examples
#' \dontrun{
#'          TCGADownload(data,'folder')
#' }
#' @export
#' @importFrom downloader download
tcgaDownload <- function(data = NULL, path = ".") {
    dir.create(path, showWarnings = FALSE)
    root <- "https://tcga-data.nci.nih.gov"
    if (!("file" %in% colnames(data))) {
        message("Downloading folders")
        for (i in 1:nrow(data)) {
            file <- paste0(path, "/", basename(data[i, "deployLocation"]))
            message(paste0("Downloading:",
                           basename(data[i, "deployLocation"])))
            if (!file.exists(file)) {
                download(paste0(root, data[i, "deployLocation"]),
                         file)
                untar(file, exdir = path)
            }
        }
    } else {
        message("Downloading files")
        for (i in 1:nrow(data)) {
            file <- paste0(path, "/", basename(data[i, "file"]))
            message(paste0("Downloading:", basename(data[i, "file"])))
            if (!file.exists(file)) {
                download(paste0(root, gsub(".tar.gz",
                                        "",
                                        data[i, "deployLocation"]),
                                "/", data[i,"file"]), file)
            }
        }
    }
}
