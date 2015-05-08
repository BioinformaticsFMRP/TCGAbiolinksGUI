
.onAttach <- function (libname, pkgname){

  file = system.file("extdata/biomics.rda", package="biomics")
  if(file.exists(file)){
    load(file,envir = as.environment("package:biomics"))
  } else {
    biomicsEnv <- as.environment("package:biomics")
    load.biosamples(biomicsEnv)
    load.systems(biomicsEnv)
    load.platforms(biomicsEnv)
    biomicsEnv$encode.db <- load.encode()
    biomicsEnv$roadmap.db <- load.roadmap()
    load.tcga(biomicsEnv)
    save(biosample.encode,biosample.roadmap,biosample.tcga,
         encode.db,platforms,roadmap.db,systems,tcga.db,platform.table,disease.table,
         file = paste0(system.file("extdata", package="biomics"),"/biomics.rda")
    )
  }
  if (!interactive() || stats::runif(1) > 0.1) return()
  welcome.message <- paste0(
    " ==================================================\n",
    " |       _   _____   _     _   _   ____      \n",
    " |____  | | |     | | \\   / | | | |      __    \n",
    " |    | | | |     | |  \\_/  | | | |     |__    \n",
    " |____| |_| |_____| |       | |_| |____  __|   \n",
    " --------------------------------------------------\n",
    " Search, download & analyse - ENCODE, TCGA, ROADMAP\n",
    " Version:",utils::packageVersion("biOMICs"),"\n",
    " ==================================================\n",
    "Use suppressPackageStartupMessages to eliminate    \n",
    "package startup messages."
  )
  packageStartupMessage(welcome.message)

}
