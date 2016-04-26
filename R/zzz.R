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
    env <- as.environment("package:TCGAbiolinksGUI")
    load.maf(env)
    maf.files <- get("maf.files", envir = env)
    use_data(maf.files, internal = TRUE, overwrite = T)
}
