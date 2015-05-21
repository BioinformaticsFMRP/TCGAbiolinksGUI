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
        load.encode(env)
        load.roadmap(env)
        load.tcga(env)
        tcga.db <- get("tcga.db")
        encode.db <- get("encode.db")
        roadmap.db <- get("roadmap.db")
        platform.table <- get("platform.table")
        disease.table <- get("disease.table")
        biosample.encode <- get("biosample.encode")
        biosample.roadmap <- get("biosample.roadmap")
        biosample.tcga <- get("biosample.tcga")
        systems <- get("systems")
        platforms <- get("platforms")

        env$success <- FALSE
        env$solution <- ""
        save(biosample.encode, biosample.roadmap, biosample.tcga,
            encode.db, platforms, roadmap.db, systems,
            tcga.db, platform.table, disease.table,
            file = file.path(system.file("extdata", package = "biOmics"),
                            "biomics.rda")
            )
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
    old <-  finf[as.numeric(difftime(Sys.time(), finf[, "mtime"],
                                     units = "days"))[1] > days,]
    return(nrow(old) == 1)
}
