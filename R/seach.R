
# internal: used by biOMICs.search
systemSearch <- function(term){
  #print (term)
  if(!success){
    found <- intersect(term, systems$BTO)
    if(length(found) > 0){
      success <<- T
      solution <<- found
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


  # search in search cache
  if(exists("search.cache")){
    idx <- grep(term, search.cache[,1])
    if(length(idx)> 0){
      message("found in cache")
      res <- search.cache[idx[1],2]
      success <<- T
      solution <<- res
      return ()
    }
  }
  # search in roamap
  # Could be - Accession, Sample.Name, Experiment
  idx <- grep(term, biosample.roadmap$biosample)
  if(length(idx)> 0){
    message("found in roadmap")
    res <- biosample.roadmap[idx[1],]$BTO
    if(!is.na(res)){
      success <<- T
      solution <<- c(res)
      return ()
    }
  }
  # search in encode
  idx <- grep(term, biosample.encode$biosample)

  if(length(idx)> 0){
    message("found in encode")
    res <- biosample.encode[idx[1],]$BTO
    if(!is.na(res)){
      success <<- T
      solution <<- res
      return ()
    }

  }
  # search in tcga

  idx <- grep(term, biosample.tcga$biosample)
  if(length(idx) > 0){
    res <- biosample.tcga[idx[1],]$BTO
    success <<- T
    solution <<- res
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
  message(paste("biOMICS is searching for:", term, "\nSearching..."))
  start.time <- Sys.time()
  success <<- F
  solution <<- NA
  exper <<- experiment

  # Step 0: verify if experiment is valid.
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
    ont <<- 'BTO'
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
      search.cache <<- rbind(c(term,solution[1]))
    } else {
      search.cache <<- rbind(search.cache,c(term,solution[1]))
    }
  }

  return(show.results())
}

# show the results to the user
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
    idx <- apply(sapply(platforms[platforms$Standard == exper, 2], function(x){grepl(pattern = x, enc.result$assay)}),1,any)
    enc.result <- enc.result[idx,]

    idx <- apply(sapply(platforms[platforms$Standard == exper, 2], function(x){grepl(pattern = x, rmap.result$Experiment)}),1,any)
    rmap.result <- rmap.result[idx,]

    idx <- apply(sapply(platforms[platforms$Standard == exper, 2], function(x){grepl(pattern = x, tcga.result$baseName)}),1,any)

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
  if(nrow(tcga.result) > 0) {
    n.barcode <- sum(sapply(tcga.result$deployStatus,function(x){length(unlist(str_split(x,",")))}))
  } else {
    n.barcode <- 0
  }
  message(paste0("|TCGA Archives: " , nrow(tcga.result)))
  message(paste0("|     Barcodes: " , n.barcode))
  message(paste0("|ENCODE : " , nrow(enc.result)))
  message(paste0("|ROADMAP: " , nrow(rmap.result)))
  message("====================================================")
  message("Some summary images were saved in searchSummary folder")

  # Preparing the output table
  colnames(rmap.result)[1:3] <- c("ID","Sample","Experiment")
  colnames(enc.result) [1:3] <- c("ID","Sample","Experiment")
  colnames(tcga.result) [c(5,11,10)] <- c("ID","Sample","Experiment")

  results <- rbind(enc.result[1:3],
                   rmap.result[1:3],
                   tcga.result[c(5,11,10)]
  )
  results$database <- c(rep("encode",  nrow(enc.result)),
                        rep("roadmap", nrow(rmap.result)),
                        rep("tcga",    nrow(tcga.result))
  )

  dir.create("searchSummary", showWarnings = F)
  # % Experiments per database
  g <- ggplot(results, aes(factor(database), fill = Experiment)) + geom_bar(position = "fill")
  ggsave(g, file="searchSummary/experiments.pdf")

  # % Samples per database
  g <- ggplot(results, aes(factor(database), fill = Sample)) + geom_bar(position = "fill")
  ggsave(g, file="searchSummary/samples.pdf")

  return(results)

}

is.experiment <- function(experiment){
  v <- c(unique(platforms$Standard), 'all')
  if(is.element(experiment, v)){
    return (TRUE)
  }
  else{
    message(paste0('ERROR: ', experiment, ' is not an experiment.\nUse:'))
    print(unique(platforms$Standard))
    return(FALSE)
  }
}

# Using the name create two collumns Platform and Disease
tcga.db.addCol <- function(x){
  tcga.db$Platform <- ""
  tcga.db$Disease <- ""
  diseases <- sapply(strsplit(biosample.tcga$biosample,split = " - "),function(x){x[1]})
  for(i in seq_along(diseases)){
    idx <- grep(diseases[i],tcga.db$baseName)
    tcga.db[idx,]$Disease <- diseases[i]
  }

  for(i in seq_along(platform.table$name)){
    idx <- grep(platform.table[i,]$name,tcga.db$baseName)
    if(length(idx)>0){
      tcga.db[idx,]$Platform <- platform.table[i,]$name
    }
  }
}
