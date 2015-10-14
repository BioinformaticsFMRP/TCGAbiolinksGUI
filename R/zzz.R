.onAttach <- function(libname, pkgname) {

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

#' Update database
#' This function will update the biomics database
#' @keywords internal
#' @importFrom devtools use_data
update <- function(){
    env <- as.environment("package:biOmics")
    load.biosamples(env)
    load.systems(env)
    load.platforms(env)
    load.encode(env)
    load.roadmap(env)
    encode.db <- get("encode.db", envir = env)
    roadmap.db <- get("roadmap.db", envir = env)
    biosample.encode <- get("biosample.encode", envir = env)
    biosample.roadmap <- get("biosample.roadmap", envir = env)
    biosample.tcga <- get("biosample.tcga", envir = env)
    systems <- get("systems", envir = env)
    platforms <- get("platforms", envir = env)

    use_data(biosample.encode,
             biosample.roadmap,
             biosample.tcga,
             encode.db,
             platforms,
             roadmap.db,
             systems,
             internal = TRUE, overwrite = T
    )

}
