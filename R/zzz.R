
.onAttach <- function (libname, pkgname){

  file = system.file("extdata/biomics.rda", package="biOMICs")
  if(file.exists(file)){
    load(file,envir = as.environment("package:biOMICs"))
  } else {
    biomicsEnv <- as.environment("package:biOMICs")
    load.biosamples(biomicsEnv)
    load.systems(biomicsEnv)
    load.platforms(biomicsEnv)
    biomicsEnv$encode.db <- load.encode()
    biomicsEnv$roadmap.db <- load.roadmap()
    load.tcga(biomicsEnv)
    save(biosample.encode,biosample.roadmap,biosample.tcga,
         encode.db,platforms,roadmap.db,systems,tcga.db,platform.table,disease.table,
         file = paste0(system.file("extdata", package="biOMICs"),"/biomics.rda")
    )
  }

  file = system.file("extdata/GRCh.rda", package="biOMICs")
  if(file.exists(file)){
    load(file,envir = as.environment("package:biOMICs"))
  }



  welcome.message <- paste0(
    " ==================================================\n",
    " |       _   _____   _     _   _   ____      \n",
    " |____  | | |     | | \\   / | | | |      __    \n",
    " |    | | | |     | |  \\_/  | | | |     |__    \n",
    " |____| |_| |_____| |       | |_| |____  __|   \n",
    " --------------------------------------------------\n",
    " Search, download & analyse - ENCODE, TCGA, ROADMAP\n",
    " Version:0.01 \n",
    " ==================================================\n"
  )
  cat("\014")
  packageStartupMessage(welcome.message)

}
