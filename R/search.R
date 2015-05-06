
# internal: used by biOMICs.search
systemSearch <- function(term){
  #print (term)
  env <- as.environment("package:biOMICs")
  if(!success){
    found <- intersect(term, systems$BTO)
    if(length(found) > 0){
      env$success <- T
      env$solution <- found
      print(solution)
      return()
    }else{
      parentTerm <- rols::parents(term[1], ont)
      parentTerm <- names(parentTerm)
      sapply(parentTerm,function(x) {systemSearch(x)})
      return ()
    }
  }
}

map.to.bto <- function(term){
  # Does the term exists in BTO?
  aux <- term
  names(aux) <- NULL
  query <- rols::olsQuery(aux,"BTO", exact = F)

  if(length(query)>0){
    # term found in BTO
    term.found <- names(query)
    sapply(term.found,
           function (x) {
             if(!success){
               systemSearch(x)
             }
           })
  }
  if(success){return()}
  else{
    ontology <- unlist(strsplit(names(term),":"))[1]
    parentTerm <- rols::parents(names(term),ontology)

    for(i in seq(parentTerm)){
      map.to.bto(parentTerm[i])
    }
  }
}


is.mapped <- function(term){
  env <- as.environment("package:biOMICs")

  # search in search cache
  if(exists("search.cache")){
    idx <- grep(term, search.cache[,1], ignore.case = T)
    if(length(idx)> 0){
      message("found in cache")
      res <- search.cache[idx[1],2]
      env$success <- T
      env$solution <- res
      return ()
    }
  }
  # search in roamap
  # Could be - Accession, Sample.Name, Experiment
  idx <- grep(term, biosample.roadmap$biosample, ignore.case = T)
  if(length(idx)> 0){
    message("found in roadmap")
    res <- biosample.roadmap[idx[1],]$BTO
    if(!is.na(res)){
      env$success <- T
      env$solution <- c(res)
      return ()
    }
  }
  # search in encode
  idx <- grep(term, biosample.encode$biosample, ignore.case = T)

  if(length(idx)> 0){
    message("found in encode")
    res <- biosample.encode[idx[1],]$BTO
    if(!is.na(res)){
      env$success <- T
      env$solution <- res
      return ()
    }

  }
  # search in tcga

  idx <- grep(term, biosample.tcga$biosample, ignore.case = T)
  if(length(idx) > 0){
    res <- biosample.tcga[idx[1],]$BTO
    env$success <- T
    env$solution <- res
    return ()
  }
}

#' @title Searches a term in TCGA, ENCODE, ROADMAP
#' @description
#'  Search a term using an ontology in the TCGA, ENCODE, ROADMAP databases
#' @param term Term to be searched
#' @import rols
#' @export
biOMICs.search  <- function(term, experiment = 'all'){
  env <- as.environment("package:biOMICs")
  message(paste("biOMICS is searching for:", term, "\nSearching..."))
  start.time <- Sys.time()
  env$success <- F
  env$solution <- NA
  env$exper <- experiment

  # Step 0: verify if term is valid.
  if(!is.valid.term(term)){
    return()
  }

  # Step 0.5: verify if experiment is valid.
  if(!is.experiment(exper)){
    return()
  }

  # Step 1: verify if term search has been mapped by us.
  if(!success){
    is.mapped(term)
  }

  # Step 2: search for the term in the BTO ontology.
  if(!success){
    message("Not found in the cache, searching in the ontology...")
    # Term has not been mapped before.
    env$ont <- 'BTO'
    query <- rols::olsQuery(term, ont, exact = T)

    if(length(query) > 0){
      term.found <- names(query)
      systemSearch(term.found)
    } else{
      query <- rols::olsQuery(term, ont, exact = F)
      if(length(query)>0){
        term.found <- names(query)
        sapply(term.found,
               function (x) {
                 if(!rols::isIdObsolete(x,ont)){
                   systemSearch(x)
                 }
               })
      }
    }

    # Not found in bto or cache
    # we are going to search in other ontologies
    # and see if from the other ontologies we can reach BTO
    if(!success){
      message("Not found in the cache, not found in BTO...")

      ont.vec <- c("EFO","CL","UBERON")
      for(i in seq(ont.vec)){
        query <- rols::olsQuery(term, ont.vec[i])
        if(length(query) > 0){
          # term was found in other ontology!
          for(i in seq(query)){
            map.to.bto(query[i])
            if(success){break}
          }
        }}
    }
  }

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message(paste("Time taken: ",round(time.taken,2),"s"))

  if(success){
    if(!exists("search.cache")){
      env$search.cache <- rbind(c(term,solution[1]))
    } else {
      env$search.cache <- rbind(search.cache,c(term,solution[1]))
    }
  }

  return(show.results())
}

# show the results to the user
#' @import ggplot2
show.results <- function(){
  # Get the samples that matches the result of the query
  # Databases were matched manually to systems
  #print(solution)

  pat <- unlist(strsplit(solution[1],","))
  idx <- apply(sapply(pat, function(x){grepl(x,biosample.encode$BTO)}),1,any)
  enc.samples  <- biosample.encode[idx,]$biosample

  idx <- apply(sapply(pat, function(x){grepl(x,biosample.roadmap$BTO)}),1,any)
  rmap.samples  <- biosample.roadmap[idx,]$biosample

  idx <- apply(sapply(pat, function(x){grepl(x,biosample.tcga$BTO)}),1,any)
  tcga.samples   <- biosample.tcga[idx,]$biosample


  # Select the database rows
  enc.result   <- encode.db [is.element(encode.db$biosample   ,enc.samples ),]
  rmap.result  <- roadmap.db[is.element(roadmap.db$Sample.Name,rmap.samples),]
  disease <- sapply(strsplit(tcga.samples,split = " - "),function(x){x[1]})
  tcga.result  <- tcga.db   [is.element(tcga.db$Disease,disease),]
  # Select experiments
  if(!exper == 'all'){
    message("Filtering by experiment")
    idx <- apply(sapply(platforms[grep(exper, platforms$Standard, ignore.case = T), 2],
                        function(x){grepl(x, enc.result$assay, ignore.case = T )}),1,any)
    enc.result <- enc.result[idx,]

    idx <- apply(sapply(platforms[grep(exper, platforms$Standard, ignore.case = T), 2],
                        function(x){grepl(x, rmap.result$Experiment, ignore.case = T)}),1,any)
    rmap.result <- rmap.result[idx,]

    idx <- apply(sapply(platforms[grep(exper, platforms$Standard, ignore.case = T), 2],
                        function(x){grepl(x, tcga.result$baseName, ignore.case = T)}),1,any)
    if(length(idx)>0){
      tcga.result <- tcga.result[idx,]
    } else {
      tcga.result <-  tcga.result[1:nrow(tcga.result),]
    }

  }

  message("============ Summary of results found ==============")
  sapply(pat, function(x){
    message(paste("|Mapped to:", subset(systems,BTO==x)$system))
  })
  message("---------- Number of terms in the system -----------")
  message(paste0("|TCGA   : " , length(tcga.samples)))
  message(paste0("|ENCODE : " , length(enc.samples)))
  message(paste0("|ROADMAP: " , length(rmap.samples)))

  message("--------------- Number of samples ------------------")
  message(paste0("|TCGA Archives: " , nrow(tcga.result)))
  message(paste0("|ENCODE : " , nrow(enc.result)))
  message(paste0("|ROADMAP: " , nrow(rmap.result)))
  message("====================================================")
  message("Some summary images were saved in searchSummary folder")

  # Preparing the output table
  colnames(rmap.result)[1:3] <- c("ID","Sample","Experiment")
  colnames(enc.result) [1:3] <- c("ID","Sample","Experiment")
  colnames(tcga.result) [c(8,12,11)] <- c("ID","Sample","Experiment")

  results <- rbind(enc.result[1:3],
                   rmap.result[1:3],
                   tcga.result[c(8,12,11)]
  )
  database <- c(rep("encode",  nrow(enc.result)),
                rep("roadmap", nrow(rmap.result)),
                rep("tcga",    nrow(tcga.result))
  )
  results <- cbind(database,results)
  dir.create("searchSummary", showWarnings = F)
  # % Experiments per database
  g <- ggplot(results, aes(factor(database), fill = Experiment)) + geom_bar(position = "fill")
  ggsave(g, file="searchSummary/experiments.pdf", height=14,scale=1.5)

  # % Samples per database
  g <- ggplot(results, aes(factor(database), fill = Sample)) + geom_bar(position = "fill")
  ggsave(g, file="searchSummary/samples.pdf", height=14,scale=1.5)

  return(results)

}

is.experiment <- function(experiment){
  v <- c(unique(platforms$Standard), 'all')
  if((length(grep(experiment, v, ignore.case = T)) > 0) & (nchar(experiment) > 3)){
    return (TRUE)
  }
  else{
    message(paste0('ERROR: ', experiment, ' is not an experiment or has less than 3 characters.\nUse:'))
    print(unique(platforms$Standard))
    return (FALSE)
  }
}

is.valid.term <- function(term){
  if(nchar(term) >= 3){
    return(TRUE)
  }
  else{
    message(paste0('ERROR: ', term, ' is not valid. Specify a term of at least 3 characters'))
    return(FALSE)
  }
}
