#' Create query from UI to be used with esearch - not completed
#' Create query from UI to be used with esearch - not completed
#' @param List of parameters parameter = (Term,Field,Connector)
#'        First parameter should contain the terms,
#'        Field NULL which means all fields and connector NULL
#'        For roadmap one parameter should be
#'        c('roadmap epigenomics','Project','AND'),
#'        For samples: c ('gsm','Filter','AND') )
#' @examples
#'  eGeo(
#'    c ('h1 cell RRBS',NULL,NULL),
#'    c ('roadmap epigenomics','Project','AND'),
#'    c ('gsm','Filter','AND') )
#'
#' @keywords internal
#' @export
#' @return Return the list to search in geo database. Used by
#'         the UI
eGeo <- function(...) {
    input_list <- list(...)
    output_list <- lapply(X = input_list, function(x) {
        if (length(x) > 2) {
            paste(x[[3]], paste(paste(x[[1]], x[[2]], sep = "["),
                                NULL, sep = "]"))
        } else {
            paste(paste(NULL, x[[1]], sep = "("), NULL, sep = ")")
        }
    })
    output_list <- paste(output_list, collapse = " ")
    return(output_list)
}

#' @title Download data from GEO database
#' @description Download data from GEO database using reutils library
#'              Input should be a GEO query, you can you this site
#'              to create the query  http://www.ncbi.nlm.nih.gov/gds/advanced
#' @param query - GEO query -  http://www.ncbi.nlm.nih.gov/gds/advanced
#' @examples
#'    geoDownloader ("(h1[All Fields] AND RRBS[All Fields])")
#'  query <- "((GSM956006) AND gsm[Filter]) AND roadmap epigenomics[Project]"
#'  query <- "((h1 cell ) AND gsm[Filter]) AND roadmap epigenomics[Project]"
#' @seealso \url{http://www.ncbi.nlm.nih.gov/geo/info/download.html}
#' @seealso \url{http://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html}
#' @seealso \url{http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch}
#' @name geoDownloader
#' @export
#' @importFrom reutils esummary content esearch
#' @importFrom XML xmlToList
#' @return Download data
geoDownloader <- function(query) {

    # search on geo databse(gds) for iQuery
    pmids <- esearch(query, "gds", retmax = 10000)

    if (pmids@.xData$no_errors()) {
        nbFiles <- as.numeric(xmlToList(pmids@.xData$get_content())$Count)

        if (nbFiles > 0) {
            # Found more than one result?

            aux <- esummary(pmids)
            info <- content(aux, "parsed")

            if (nbFiles > 1) {
                ftps <- lapply(info, function(x) x$FTPLink)
            } else {
                ftps <- info$FTPLink
            }
            geoDownloaderLinks(ftps)
        }
    }
}

#' @title Verify user input
#' @description Verify user input existis in roadmap
#' @param input - user input of samples
#' @param list - roadmap list of experiments or samples
#' @examples
#' \dontrun{
#'  validateInput(samples,roadmap$samplesList)
#' }
#' @keywords internal
#' @return If uinput is valid or not
validateInput <- function(input = NULL, list = NULL) {
    for (i in seq_along(input)) {
        if ( !(input[i] %in% list)) {
            print(paste0("ERROR: ", input[i], " does not exists in roadmap"))
            print(paste0("Possible values can be seen in objects: ",
                        "roadmap$experimentList and roadmap$samplesLis "))
            # print(list)
        }
    }
}

#' @title Filter roadmap data
#' @description Select data from roadmap
#' @param samples - list of samples
#' @param experiments - list of experiments
#' @examples
#'    link <- filterRmapData (
#'     samples = c('H1 cell line','H9 cell line'),
#'     experiments = c('RRBS','H4K91ac')
#'    )
#' @return List of ftp to download data
#' @export
#' @return Get list of files doing the samples x experiments
filterRmapData <- function(samples = NULL, experiments = NULL) {
    # Load ROAMP table if it doesn't exists
    roadmap.db <- getOption("roadmap.db")

    if (!is.null(samples))
        validateInput(samples, roadmap.db$Sample.Name)
    if (!is.null(experiments))
        validateInput(experiments, roadmap.db$Experiment)
    # get list of ftps use outer in order to combine each element
    # of x and y
    if (!is.null(samples) && !is.null(experiments)) {
        link <- unlist(outer(X = samples,
                            Y = experiments,
                            Vectorize(function(x,y) {
                            as.character(roadmap.db$GEO.FTP[
                                roadmap.db$Sample.Name == x &
                                roadmap.db$Experiment == y])
                             })))
    } else if (!is.null(samples)) {
        link <- unlist(lapply(samples, function(x) {
            as.character(roadmap.db$GEO.FTP[roadmap.db$Sample.Name == x])
        }))
    } else if (!is.null(experiments)) {
        link <- unlist(lapply(experiments, function(x) {
            as.character(roadmap.db$GEO.FTP[roadmap.db$Sample.Name == x])
        }))
    }
    invisible(link)
}

#' Download all files of ftp directory
#' Download all files of ftp directory
#' @param iFTP: ftp directory adress
#' @param  iFileName - list of files in the ftp path
#'        (use .getFileNames to get the list)
#' @param iOut - folder where files will be downloaded
#' @return compressedFiles - compressed files downloaded
#' @examples
#'  downloadFromGEO (
#'   'ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM409nnn/GSM409307/suppl/',
#'    c ('GSM409307_UCSD.H1.H3K4me1.LL228.bed.gz'),
#'    'path_to_folder_GSM409307' )
#' @keywords internal
#' @importFrom pbapply pblapply
#' @importFrom downloader download
downloadFromGEO <- function(iFTP, iFileName, iOut, gui = FALSE) {
    compresssedFiles <- c()
    if (gui)
        setOptionsProgressBar(title = "Roadmap data", label = "Downloading")
    gzFiles <- pblapply(iFileName, function(x) {
        if (!file.exists(paste0(iOut, "/", x)))
            download(paste0(iFTP, x), paste0(iOut, "/", x))
        compresssedFiles <- c(compresssedFiles, file)
    })
    return(gzFiles)
}

#' Download a files from a list o ftp directory and uncompress it
#' Download a files from a list o ftp directory and uncompress it
#' @param iFTPs: list of ftp directory adress
#'        iOut - folder where files will be downloaded
#' @examples
#'  geoDownloaderLinks (
#'     c('ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM409nnn/GSM409307/suppl/',
#'       'ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM409nnn/GSM409312/suppl/'
#'       ),
#'     'path_to_folder_GSM409307'
#'   )
#' @keywords internal
#' @export
#' @return download geo links
geoDownloaderLinks <- function(iFTPs, iOut, gui = FALSE) {
    nbFiles <- length(iFTPs)
    if (nbFiles > 0) {
        # Found more than one result?
        dir.create(iOut)

        for (ftp in iFTPs) {
            if (ftp != "") {
                # create folder for experiment
                dirName <- tail(unlist(strsplit(ftp, "/")), n = 2)[1]
                outPath <- paste(iOut, dirName, sep = "/")
                dir.create(outPath)
                # get list of files names in the ftp
                fileName <- getFileNames(ftp)
                compresssedFiles <- downloadFromGEO(ftp, fileName,
                                                    outPath, gui)
                print("Uncompressing")
                uncompress(compresssedFiles)

            }
        }
        # prepareInfoTable (ret$link,iOut)
    }
}

# Uncompress a list of files - it does not remove the
# compressed file @keywords internal
#' @importFrom R.utils gunzip
uncompress <- function(iFiles) {
    if (!is.null(iFiles)) {
        pbapply::pblapply(iFiles, function(x) {
            if (file.exists(x) & tools::file_ext(x) == "gz") {
                gunzip(x, remove = FALSE, skip = TRUE,
                        ext = "gz")
                print(paste0("Uncompressed: ", x))
            }
        })
    } else {
        return(NULL)
    }
}

# Show ftp links downloaded @keywords internal
prepareInfoTable <- function(iLink, iOut) {
    if (length(iLink) > 0) {
        df <- do.call(rbind.data.frame, iLink)
        result$g_nbFiles <- nrow(df)
        names(df) <- paste0("Files downloaded into:", getwd(),
                            "/", iOut)
        result$df <- df
    }
}

# Get all files in the ftp directory @keywords internal
getFileNames <- function(ftp) {
    filePath <- unlist(strsplit(RCurl::getURL(ftp, ftp.use.epsv = FALSE,
                                            dirlistonly = TRUE), " "))
    # remove \r\n from name
    fileName <- as.list(strsplit(tail(filePath, n = 1), "\r*\n")[[1]])
}

#' @title Verify user input
#' @description Verify user input existis in roadmap
#' @param input - user input of samples
#' @param list - roadmap list of experiments or samples
#' @examples
#' \dontrun{
#'  validateInput(samples,roadmap$samplesList)
#' }
#' @keywords internal
validateInput <- function(input = NULL, list = NULL) {
    for (i in seq_along(input)) {
        if ( !(input[i] %in% list)) {
            print(paste0("ERROR: ", input[i], " does not exists in roadmap"))
            print(paste0("Possible values can be seen in objects: ",
                        "roadmap$experimentList and roadmap$samplesLis "))
            # print(list)
        }
    }
}

#' @title Filter roadmap data
#' @description Select data from roadmap
#' @param samples - list of samples
#' @param experiments - list of experiments
#' @examples
#' \dontrun{
#'    link <- filterRmapData (
#'     samples = c('H1 cell line','H9 cell line'),
#'     experiments = c('RRBS','H4K91ac')
#'    )
#' }
#' @return List of ftp to download data
#' @export
filterRmapData <- function(samples = NULL, experiments = NULL) {
    # Load ROAMP table if it doesn't exists
    roadmap.db <- getOption("roadmap.db")

    if (!is.null(samples))
        validateInput(samples, roadmap.db$Sample.Name)
    if (!is.null(experiments))
        validateInput(experiments, roadmap.db$Experiment)
    # get list of ftps use outer in order to combine each element
    # of x and y
    if (!is.null(samples) && !is.null(experiments)) {
        link <- unlist(outer(X = samples,
                            Y = experiments,
                            Vectorize(function(x,y) {
                            as.character(roadmap.db$GEO.FTP[
                                roadmap.db$Sample.Name == x &
                                roadmap.db$Experiment == y])
                             })))
    } else if (!is.null(samples)) {
        link <- unlist(lapply(samples, function(x) {
            as.character(roadmap.db$GEO.FTP[roadmap.db$Sample.Name == x])
        }))
    } else if (!is.null(experiments)) {
        link <- unlist(lapply(experiments, function(x) {
            as.character(roadmap.db$GEO.FTP[roadmap.db$Sample.Name == x])
        }))
    }
    invisible(link)
}

#' Download all files of ftp directory
#' Download all files of ftp directory
#' @param iFTP: ftp directory adress
#' @param  iFileName - list of files in the ftp path
#'        (use .getFileNames to get the list)
#' @param iOut - folder where files will be downloaded
#' @return compressedFiles - compressed files downloaded
#' @examples
#' \dontrun{
#'  downloadFromGEO (
#'   'ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM409nnn/GSM409307/suppl/',
#'    c ('GSM409307_UCSD.H1.H3K4me1.LL228.bed.gz'),
#'    'path_to_folder_GSM409307' )
#' }
#' @keywords internal
#' @importFrom pbapply pblapply
#' @importFrom downloader download
downloadFromGEO <- function(iFTP, iFileName, iOut, gui = FALSE) {
    compresssedFiles <- c()
    if (gui)
        setOptionsProgressBar(title = "Roadmap data", label = "Downloading")
    gzFiles <- pblapply(iFileName, function(x) {
        if (!file.exists(paste0(iOut, "/", x)))
            download(paste0(iFTP, x), paste0(iOut, "/", x))
        compresssedFiles <- c(compresssedFiles, file)
    })
    return(gzFiles)
}

#' Download a files from a list o ftp directory and uncompress it
#' Download a files from a list o ftp directory and uncompress it
#' @param iFTPs: list of ftp directory adress
#'        iOut - folder where files will be downloaded
#' @examples
#' \dontrun{
#'  geoDownloaderLinks (
#'     c('ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM409nnn/GSM409307/suppl/',
#'       'ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM409nnn/GSM409312/suppl/'
#'       ),
#'     'path_to_folder_GSM409307'
#'   )
#' }
#' @keywords internal
#' @export
geoDownloaderLinks <- function(iFTPs, iOut, gui = FALSE) {
    nbFiles <- length(iFTPs)
    if (nbFiles > 0) {
        # Found more than one result?
        dir.create(iOut)

        for (ftp in iFTPs) {
            if (ftp != "") {
                # create folder for experiment
                dirName <- tail(unlist(strsplit(ftp, "/")), n = 2)[1]
                outPath <- paste(iOut, dirName, sep = "/")
                dir.create(outPath)
                # get list of files names in the ftp
                fileName <- getFileNames(ftp)
                compresssedFiles <- downloadFromGEO(ftp, fileName,
                                                    outPath, gui)
                print("Uncompressing")
                uncompress(compresssedFiles)

            }
        }
        # prepareInfoTable (ret$link,iOut)
    }
}

# Uncompress a list of files - it does not remove the
# compressed file @keywords internal
#' @importFrom R.utils gunzip
uncompress <- function(iFiles) {
    if (!is.null(iFiles)) {
        pbapply::pblapply(iFiles, function(x) {
            if (file.exists(x) & tools::file_ext(x) == "gz") {
                gunzip(x, remove = FALSE, skip = TRUE,
                        ext = "gz")
                print(paste0("Uncompressed: ", x))
            }
        })
    } else {
        return(NULL)
    }
}

# Show ftp links downloaded @keywords internal
prepareInfoTable <- function(iLink, iOut) {
    if (length(iLink) > 0) {
        df <- do.call(rbind.data.frame, iLink)
        result$g_nbFiles <- nrow(df)
        names(df) <- paste0("Files downloaded into:", getwd(),
                            "/", iOut)
        result$df <- df
    }
}

# Get all files in the ftp directory @keywords internal
getFileNames <- function(ftp) {
    filePath <- unlist(strsplit(RCurl::getURL(ftp, ftp.use.epsv = FALSE,
                                            dirlistonly = TRUE), " "))
    # remove \r\n from name
    fileName <- as.list(strsplit(tail(filePath, n = 1), "\r*\n")[[1]])
}
