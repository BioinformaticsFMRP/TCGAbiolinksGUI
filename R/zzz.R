
.onAttach <- function (libname, pkgname){
  options(search.cache=NULL)
  file = system.file("extdata/biomics.rda", package="biOmics")
  if(file.exists(file)){
    load(file,envir = as.environment("package:biOmics"))
  } else {
    biomicsEnv <- as.environment("package:biOmics")
    load.biosamples(biomicsEnv)
    load.systems(biomicsEnv)
    load.platforms(biomicsEnv)
    biomicsEnv$encode.db <- load.encode()
    biomicsEnv$roadmap.db <- load.roadmap()
    load.tcga(biomicsEnv)
    biomicsEnv$success <- F
    biomicsEnv$solution <- ""
    save(biosample.encode,biosample.roadmap,biosample.tcga,
         encode.db,platforms,roadmap.db,systems,tcga.db,platform.table,disease.table,
         file = paste0(system.file("extdata", package="biOmics"),"/biomics.rda")
    )
  }


  file = system.file("extdata/GRCh.rda", package="biOmics")
  if(file.exists(file)){
    load(file,envir = as.environment("package:biOmics"))
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
    " Version:",utils::packageVersion("biOmics"),"\n",
    " ==================================================\n",
    " Use suppressPackageStartupMessages to eliminate    \n",
    " package startup messages."
  )
  packageStartupMessage(welcome.message)

}
