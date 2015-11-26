#' @title biOmicsDownload
#' @description Download data selected using the biOmics.search function
#' @param lines biOmicsSearch output
#' @param enc.file.type Extension to be downloaed from encode database
#' @param rmap.file.type Extension to be downloaed from roadmap database
#' @seealso biOmicsSearch
#' @export
#' @return Saves the file of the line files
#' @examples
#' query <- biOmicsSearch("u87")
#' biOmicsDownload(query, enc.file.type = "bam",rmap.file.type = "bed")
biOmicsDownload <- function(lines=NULL,
                            enc.file.type = NULL,
                            rmap.file.type = NULL) {

    if (is.null(lines)) stop("Please set lines parameter")

    encode.lines <- subset(lines, lines$database == "encode")
    rmap.lines <- subset(lines, lines$database == "roadmap")
    tcga.lines <- subset(lines, lines$database == "tcga")

    # -------------- download Encode
    if (dim(encode.lines)[1] > 0) {
        message("==== Encode download ====")
        encode.lines <- subset(encode.db,
                               encode.db$accession == encode.lines$ID)
        #encodeDownload(encode.lines, enc.file.type)
    }

    # -------------- download ROADMAP
    if (dim(rmap.lines)[1] > 0) {
        message("==== Roadmap download ====")
        rmap.lines <- subset(roadmap.db,
                             roadmap.db$EID == rmap.lines$ID)
        roadmapDownload(rmap.lines, rmap.file.type)
    }
    # ---------------- download TCGA TODO: add filters, folder to

    # This will be replaced by TCGAbiolinks
    if (dim(tcga.lines)[1] > 0) {
        message("==== TCGA download ====")
        tcga.db <- TCGAquery()
        tcga.lines <- tcga.db[tcga.db$name == tcga.lines$ID,]
        TCGAdownload(tcga.lines, path = "TCGA")
    }


}

#' @title Download encode data
#' @description Download encode data selected using the encodeSearch
#' @param lines encode.search output
#' @param path Folder to save the file
#' @param type extesion of files to be downloaded
#' @export
#' @examples
#'  query <- encodeSearch(biosample = "brain")
#'  encodeDownload(query,type = "bam",path = "encode")
#' @importFrom downloader download
#' @importFrom RCurl getURL
#' @importFrom rjson fromJSON
#' @return Download encode data into path
encodeDownload <- function(lines, type = NULL, path = ".") {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)

    encode.url <- "https://www.encodeproject.org/"
    json <- "/?format=JSON"

    for (i in 1:dim(lines)[1]) {
        id <- lines$accession[i]
        url <- paste0(encode.url, "experiments/", id, json)
        dir.create(file.path(path,id), showWarnings = FALSE, recursive = TRUE)
        # get list of files
        item <- fromJSON(getURL(url, dirlistonly = TRUE,
                                .opts = list(ssl.verifypeer = FALSE)
        )
        )[["files"]]

        files <- sapply(item, function(x) { x$href })

        if (!is.null(type)) {
            idx <- unlist(lapply(type, function(x) {
                grep(x, files)
            }))
            files <- files[idx]
        }
        # download files
        for (j in files) {
            link <- gsub("https:","http:",paste0(encode.url, j))
            fileout <- file.path(path, id, basename(link))

            # Downloader library is not working here =/
            if(!file.exists(fileout))
                download.file(link, fileout, method = "wget")
        }
    }

}

#' @title Download roadmap data
#' @description Download roadmap data selected using the roadmap.search
#' @param lines roadmap.search output
#' @param path Folder to save the file
#' @param type extesion of files to be downloaded
#' @export
#' @import AnnotationHub
#' @importFrom stringr str_replace_all
#' @return Download romapdata into path
#' @examples
#' query <- roadmapSearch(sample = "H1 cell line", experiment = "RRBS")
#' roadmapDownload(query,type = "bed", path = "roadmap")
roadmapDownload <- function(lines, type = NULL, path = ".") {
    ah = AnnotationHub()
    for (i in 1:dim(lines)[1]) {
        epiFiles <- query(ah, c("EpigenomeRoadMap", lines[i,]$EID))
        for(j in names(epiFiles@.db_uid)){
            file <- epiFiles[j]$sourceurl
            download(file,basename(file))
        }
    }
}
