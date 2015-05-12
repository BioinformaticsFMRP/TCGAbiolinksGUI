formatTarget <- function(target) {
    if (!is.null(target)) {
        return(paste("target.investigated_as", gsub(" ", "%20",
            target), sep = "=", collapse = "&"))
    }
    return(NULL)
}

formatSample <- function(sample) {
    if (!is.null(sample)) {
        return(paste("replicates.library.biosample.biosample_type",
            gsub(" ", "%20", sample), sep = "=", collapse = "&"))
    }
    return(NULL)
}

formatFType <- function(type) {
    if (!is.null(type)) {
        return(paste("files.file_format", type, sep = "=",
            collapse = "&"))
    }
    return(NULL)
}

formatAssay <- function(assay) {
    if (!is.null(assay)) {
        return(paste("assay_term_name", assay, sep = "=", collapse = "&"))
    }
    return(NULL)
}

formatAssembly <- function(assambly) {
    if (!is.null(assambly)) {
        return(paste("assembly", assambly, sep = "=", collapse = "&"))
    }
    return(NULL)
}

formatSearch <- function(seach) {
    if (!is.null(seach))
        return(paste("searchTerm", gsub(" ", "+", seach), sep = "="))
    return(NULL)
}

formatType <- function(type) {
    if (!is.null(type))
        return(paste("type", type, sep = "="))
    return(NULL)
}

#' @title Download data from ENCODDE project
#' @description Downloads data from ENCODE project
#' @param search - bone chip
#' @param type -   'experiment' 'publication' 'software' ''antibody_lot' 'web'
#' @param target - 'transcription factor' 'RNA binding protein'
#'                  'tag' 'histone modification'
#' @param sample - 'tissue' 'primary cell'
#' @param ftype  - 'bam' 'bigWig' 'bed_broadPeak' 'broadPeak' 'fastq'
#' @param assay - 'ChIP-seq' 'RIP-chip' 'Repli-chip'
#' @param assembly - 'hg19' 'mm9'
#' @param out - path to save files
#' @examples
#' \dontrun{
#'    encodeDownloader('bone chip',
#'                    'experiment',
#'                    'transcription factor',
#'                    'primary cell',
#'                    'bigWig',
#'                    'ChIP-seq',
#'                    'mm9',
#'                    'diretory_name_to_save_files'
#'                    )
#' }
#' @seealso \url{https://www.encodeproject.org/search/}
#' @seealso \url{https://www.encodeproject.org/help/rest-api/}
#' @name encodeDownloader
#' @export
#' @return Save encode data into folder, return data downloaded
encodeDownloader <- function(search, type, target, sample,
    ftype, assay, assembly, out) {
    # Constant parameters
    encodePath <- "https://www.encodeproject.org/"
    searchPath <- "search/"
    format <- "format=json"
    frame <- "frame=embedded"

    searchTerm <- formatSearch(search)
    type <- formatType(type)
    target <- formatTarget(target)
    biosample <- formatSample(sample)
    fType <- formatFType(ftype)
    assay <- formatAssay(assay)
    assembly <- formatAssembly(assembly)

    URLoptions <- paste(searchTerm, format, frame, target, type,
        biosample, fType, assay, assembly, sep = "&")
    URL <- paste0(encodePath, paste(searchPath, URLoptions, sep = "?"))

    data <- rjson::fromJSON(RCurl::getURL(URL, dirlistonly = TRUE,
        .opts = list(ssl.verifypeer = FALSE)), unexpected.escape = "keep")

    if (data$total > 0) {
        # Output folder
        dir.create(out, showWarnings = FALSE)
        df <- data.frame()
        # Downloads all files from the search
        for (j in 1:length(data$"@graph")) {
            nbFiles <- length(data$"@graph"[[j]]$files)
            print(paste0("Downloading experiment ", j))

            for (i in 1:nbFiles) {
                filePath <- data$"@graph"[[j]]$files[[i]]$href
                dn <- paste0(encodePath, filePath)
                fileOut <- paste(out, unlist(strsplit(data$"@graph"[[j]]$"@id",
                  "/"))[3], sep = "/")
                if (!file.exists(fileOut)) {
                  dir.create(fileOut, showWarnings = FALSE)
                }
                fileOut <- paste(fileOut, unlist(strsplit(filePath,
                  "/"))[5], sep = "/")

                # Did the user select a file format
                if (is.null(type) ||
                    data$"@graph"[[j]]$files[[i]]$file_format %in% type){
                    downloader::download(dn, fileOut, quiet = TRUE)
                  print(paste0("Downloaded file: ", fileOut))

                }
                df <- rbind(df, data.frame(dn))
            }
        }
        names(df) <- paste0("Files downloaded into:", getwd(),
            "/", out)
        result$df <- df
    } else {
        print(data$notification)
    }
}
