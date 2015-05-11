.onAttach <- function(libname, pkgname) {
    options(search.cache = NULL)
    file <- system.file("extdata/biomics.rda", package = "biOmics")
    if (file.exists(file) && !is.old(file)) {
        load(file, envir = as.environment("package:biOmics"))
    } else {
        env <- as.environment("package:biOmics")
        load.biosamples(env)
        load.systems(env)
        load.platforms(env)
        env$encode.db <- load.encode()
        env$roadmap.db <- load.roadmap()
        load.tcga(env)
        env$success <- FALSE
        env$solution <- ""
        save(biosample.encode, biosample.roadmap, biosample.tcga,
            encode.db, platforms, roadmap.db, systems, tcga.db,
            platform.table, disease.table, file = paste0(system.file("extdata",
                package = "biOmics"), "/biomics.rda"))
    }

    file <- system.file("extdata/GRCh.rda", package = "biOmics")
    if (file.exists(file)) {
        load(file, envir = as.environment("package:biOmics"))
    }

    if (!interactive() || stats::runif(1) > 0.1)
        return()
    welcome.message <- paste0(
        " ==================================================\n",
        " |       _   _____   _     _   _   ____      \n",
        " |____  | | |     | | \\   / | | | |      __    \n",
        " |    | | | |     | |  \\_/  | | | |     |__    \n",
        " |____| |_| |_____| |       | |_| |____  __|   \n",
        " --------------------------------------------------\n",
        " Search, download & analyse - ENCODE, TCGA, ROADMAP\n",
        " Version:", utils::packageVersion("biOmics"), "\n",
        " ==================================================\n",
        " Use suppressPackageStartupMessages to eliminate    \n",
        " package startup messages.")
    packageStartupMessage(welcome.message)
}

is.old <- function(file = NULL, days = 20) {
    finf <- file.info(file, extra_cols = FALSE)
    return(nrow(finf[difftime(Sys.time(), finf[, "mtime"], units = "days") >
        days, 1:4]) == 1)
}
