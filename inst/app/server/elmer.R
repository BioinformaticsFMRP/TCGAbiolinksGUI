#------------------------------------------------
# ELMER
# -----------------------------------------------
#--------------------- START controlling show/hide states -----------------
#shinyjs::hide("survivalplotgroup")
#shinyjs::hide("survivalplotMain")
#shinyjs::hide("survivalplotLegend")
#shinyjs::hide("survivalplotLimit")
#shinyjs::hide("survivalplotPvalue")
observeEvent(input$scatter.plot.type, {
    type <- isolate({input$elmerPlotType})
    if(type =="scatter.plot"){
        scatter.type <- isolate({input$scatter.plot.type})
        if(scatter.type == "tf"){
            shinyjs::show("scatter.plot.tf")
            shinyjs::show("scatter.plot.motif")
            shinyjs::hide("scatter.plot.genes")
            shinyjs::hide("scatter.plot.probes")
            shinyjs::hide("scatter.plot.nb.genes")
        } else if(scatter.type == "pair"){
            shinyjs::hide("scatter.plot.tf")
            shinyjs::hide("scatter.plot.motif")
            shinyjs::show("scatter.plot.genes")
            shinyjs::show("scatter.plot.probes")
            shinyjs::hide("scatter.plot.nb.genes")
        } else {
            shinyjs::hide("scatter.plot.tf")
            shinyjs::hide("scatter.plot.motif")
            shinyjs::hide("scatter.plot.genes")
            shinyjs::show("scatter.plot.probes")
            shinyjs::show("scatter.plot.nb.genes")
        }
    }
})
observeEvent(input$elmerPlotType, {
    type <- isolate({input$elmerPlotType})
    if(type =="scatter.plot"){
        shinyjs::show("scatter.plot.type")
        shinyjs::hide("schematic.plot.type")
        shinyjs::hide("schematic.plot.genes")
        shinyjs::hide("schematic.plot.probes")
        shinyjs::hide("ranking.plot.motif")
        shinyjs::hide("elmer.heatmap.annotations")
        shinyjs::hide("motif.enrichment.plot.summary")
        shinyjs::hide("ranking.plot.tf")
        shinyjs::hide("motif.enrichment.plot.or")
        shinyjs::hide("motif.enrichment.plot.loweror")
        scatter.type <- isolate({input$scatter.plot.type})
        if(scatter.type == "tf"){
            shinyjs::show("scatter.plot.tf")
            shinyjs::show("scatter.plot.motif")
            shinyjs::hide("scatter.plot.genes")
            shinyjs::hide("scatter.plot.probes")
            shinyjs::hide("scatter.plot.nb.genes")
        } else if(scatter.type == "pair"){
            shinyjs::hide("scatter.plot.tf")
            shinyjs::hide("scatter.plot.motif")
            shinyjs::show("scatter.plot.genes")
            shinyjs::show("scatter.plot.probes")
            shinyjs::hide("scatter.plot.nb.genes")
        } else {
            shinyjs::hide("scatter.plot.tf")
            shinyjs::hide("scatter.plot.motif")
            shinyjs::hide("scatter.plot.genes")
            shinyjs::show("scatter.plot.probes")
            shinyjs::show("scatter.plot.nb.genes")
        }
    } else if(type =="schematic.plot"){
        shinyjs::hide("scatter.plot.type")
        shinyjs::hide("scatter.plot.tf")
        shinyjs::hide("scatter.plot.motif")
        shinyjs::hide("elmer.heatmap.annotations")
        shinyjs::hide("scatter.plot.genes")
        shinyjs::hide("scatter.plot.probes")
        shinyjs::hide("scatter.plot.nb.genes")
        shinyjs::show("schematic.plot.type")
        shinyjs::hide("ranking.plot.motif")
        shinyjs::hide("ranking.plot.tf")
        shinyjs::hide("motif.enrichment.plot.summary")
        shinyjs::hide("motif.enrichment.plot.or")
        shinyjs::hide("motif.enrichment.plot.loweror")
        type <- isolate({input$schematic.plot.type})
        if(type =="genes"){
            shinyjs::hide("schematic.plot.probes")
            shinyjs::show("schematic.plot.genes")
        } else {
            shinyjs::show("schematic.plot.probes")
            shinyjs::hide("schematic.plot.genes")
        }
    } else if(type =="ranking.plot"){
        shinyjs::hide("scatter.plot.type")
        shinyjs::hide("scatter.plot.tf")
        shinyjs::hide("scatter.plot.motif")
        shinyjs::hide("scatter.plot.genes")
        shinyjs::hide("elmer.heatmap.annotations")
        shinyjs::hide("scatter.plot.probes")
        shinyjs::hide("scatter.plot.nb.genes")
        shinyjs::hide("schematic.plot.type")
        shinyjs::hide("schematic.plot.genes")
        shinyjs::hide("schematic.plot.probes")
        shinyjs::show("ranking.plot.motif")
        shinyjs::show("ranking.plot.tf")
        shinyjs::hide("motif.enrichment.plot.summary")
        shinyjs::hide("motif.enrichment.plot.or")
        shinyjs::hide("motif.enrichment.plot.loweror")
    } else if(type =="motif.enrichment.plot"){
        shinyjs::hide("scatter.plot.type")
        shinyjs::hide("scatter.plot.tf")
        shinyjs::hide("scatter.plot.motif")
        shinyjs::hide("scatter.plot.genes")
        shinyjs::hide("elmer.heatmap.annotations")
        shinyjs::hide("scatter.plot.probes")
        shinyjs::hide("scatter.plot.nb.genes")
        shinyjs::hide("schematic.plot.type")
        shinyjs::hide("schematic.plot.genes")
        shinyjs::hide("schematic.plot.probes")
        shinyjs::hide("ranking.plot.motif")
        shinyjs::hide("ranking.plot.tf")
        shinyjs::show("motif.enrichment.plot.summary")
        shinyjs::show("motif.enrichment.plot.or")
        shinyjs::show("motif.enrichment.plot.loweror")
    } else if(type =="heatmap.plot"){
        shinyjs::hide("scatter.plot.type")
        shinyjs::hide("scatter.plot.tf")
        shinyjs::hide("scatter.plot.motif")
        shinyjs::hide("scatter.plot.genes")
        shinyjs::show("elmer.heatmap.annotations")
        shinyjs::hide("scatter.plot.probes")
        shinyjs::hide("motif.enrichment.plot.summary")
        shinyjs::hide("scatter.plot.nb.genes")
        shinyjs::hide("schematic.plot.type")
        shinyjs::hide("schematic.plot.genes")
        shinyjs::hide("schematic.plot.probes")
        shinyjs::hide("ranking.plot.motif")
        shinyjs::hide("ranking.plot.tf")
        shinyjs::hide("motif.enrichment.plot.or")
        shinyjs::hide("motif.enrichment.plot.loweror")
    } else if(type =="volcano.plot"){
        shinyjs::hide("scatter.plot.type")
        shinyjs::hide("scatter.plot.tf")
        shinyjs::hide("scatter.plot.motif")
        shinyjs::hide("scatter.plot.genes")
        shinyjs::hide("elmer.heatmap.annotations")
        shinyjs::hide("scatter.plot.probes")
        shinyjs::hide("motif.enrichment.plot.summary")
        shinyjs::hide("scatter.plot.nb.genes")
        shinyjs::hide("schematic.plot.type")
        shinyjs::hide("schematic.plot.genes")
        shinyjs::hide("schematic.plot.probes")
        shinyjs::hide("ranking.plot.motif")
        shinyjs::hide("ranking.plot.tf")
        shinyjs::hide("motif.enrichment.plot.or")
        shinyjs::hide("motif.enrichment.plot.loweror")

    }
})

observeEvent(input$schematic.plot.type, {
    type <- isolate({input$elmerPlotType})
    if(type =="schematic.plot"){
        type <- isolate({input$schematic.plot.type})
        if(type =="genes"){
            shinyjs::hide("schematic.plot.probes")
            shinyjs::show("schematic.plot.genes")
        } else {
            shinyjs::show("schematic.plot.probes")
            shinyjs::hide("schematic.plot.genes")
        }
    }
})
#----------------------- END controlling show/hide states -----------------
observe({
    shinyFileChoose(input, 'elmermaefile', roots=get.volumes(input$workingDir), session=session,
                    restrictions=system.file(package='base'), filetypes=c('', 'rda'))
    shinyFileChoose(input, 'elmerresultsfile', roots=get.volumes(input$workingDir), session=session,
                    restrictions=system.file(package='base'), filetypes=c('', 'rda'))
    shinyFileChoose(input, 'elmermetfile', roots=get.volumes(input$workingDir), session=session,
                    restrictions=system.file(package='base'), filetypes=c('', 'rda'))
    shinyFileChoose(input, 'elmerexpfile', roots=get.volumes(input$workingDir), session=session,
                    restrictions=system.file(package='base'), filetypes=c('', 'rda'))
    updateTextInput(session, "maesavefilename",
                    value = paste0("mae",
                                   paste(input$elmermaetype,
                                         gsub("[[:punct:]]| ","_",input$elmermaesubtype),
                                         gsub("[[:punct:]]| ","_",input$elmermaesubtype2),
                                         sep = "_"),".rda"))
})


# mae create
observe({
    data <- maedata()
    updateSelectizeInput(session, 'elmermaetype', choices = {
        if(!is.null(data)) as.character(colnames(colData(data)))
    }, server = TRUE)
})
observeEvent(input$elmermaetype, {
    data <- maedata()
    updateSelectizeInput(session, 'elmermaesubtype', choices = {
        if(!is.null(data)) as.character(unique(colData(data)[,input$elmermaetype]))
    }, server = TRUE)
})
observeEvent(input$elmermaetype, {
    data <- maedata()
    updateSelectizeInput(session, 'elmermaesubtype2', choices = {
        if(!is.null(data)) as.character(unique(colData(data)[,input$elmermaetype]))
    }, server = TRUE)
})

observe({
    input$elmermaesubtype2
    input$elmermaesubtype
    group1 <- if_else(str_length(isolate({input$elmermaesubtype})) > 0,isolate({input$elmermaesubtype}), "group 1")
    group2 <- if_else(str_length(isolate({input$elmermaesubtype2})) > 0,isolate({input$elmermaesubtype2}), "group 2")
    choices <- c("hyper","hypo")
    names(choices) <-  c(paste("Probes hypermethylated in", group1,"compared to", group2),
                         paste("Probes hypomethylated in", group1,"compared to", group2))
    updateSelectizeInput(session, 'elmerdirection', choices = {
        choices
    }, server = TRUE)
})

elmer.exp <-  reactive({
    inFile <- input$elmerexpfile
    if (is.null(inFile) || inFile == 0) return(NULL)
    file  <- as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), inFile)$datapath)

    withProgress(message = 'Loading data',
                 detail = 'This may take a while...', value = 0, {
                     exp <- get(load(file))
                     incProgress(1, detail = "Completed")
                 })

    if(!class(exp) %in% c(class(as(SummarizedExperiment(),"RangedSummarizedExperiment")),class(matrix()),class(data.frame()))){
        createAlert(session, "elmermessage", "elmerAlert", title = "Data input error", style =  "danger",
                    content = paste0("Sorry, but I'm expecting a Summarized Experiment, a data frame or a Matrix object, but I got a: ",
                                     class(exp)), append = FALSE)
        return(NULL)
    } else {
        if(!all(grepl("ENG", rownames(exp)))){
            if("ensembl_gene_id" %in% names(values(exp))) rownames(exp) <- rowRanges(exp)$ensembl_gene_id
        } else {
            createAlert(session, "elmermessage", "elmerAlert", title = "Data input error", style =  "danger",
                        content = "Sorry, but I'm expecting a Summarized Experiment, a data frame or a Matrix object with EMSEMBL GENE ID as row names",
                        append = FALSE)
            return(NULL)
        }
    }
    return(exp)
})
elmer.met <-  reactive({
    inFile <- input$elmermetfile
    if (is.null(inFile) || inFile == 0) return(NULL)
    file  <- as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), inFile)$datapath)

    withProgress(message = 'Loading data',
                 detail = 'This may take a while...', value = 0, {
                     met <- get(load(file))
                     incProgress(1, detail = "Completed")
                 })

    if(!class(met) %in% c(class(as(SummarizedExperiment(),"RangedSummarizedExperiment")),class(matrix()),class(data.frame()))){
        createAlert(session, "elmermessage", "elmerAlert", title = "Data input error", style =  "danger",
                    content = paste0("Sorry, but I'm expecting a Summarized Experiment, a data frame or a Matrix object, but I got a: ",
                                     class(met)), append = FALSE)
        return(NULL)
    }

    return(met)
})


createMae <-  reactive({
    input$elmercreatemae
    exp <- elmer.exp()
    met <- elmer.met()

    if(is.null(exp)) {
        createAlert(session, "elmerinputmessage", "elmerAlert", title = "ERROR", style =  "danger",
                    content =   "Please select gene expression object", append = TRUE)
        return(NULL)
    }
    if(is.null(met)) {
        createAlert(session, "elmerinputmessage", "elmerAlert", title = "ERROR", style =  "danger",
                    content =   "Please select DNA methylation object", append = TRUE)
        return(NULL)
    }

    withProgress(message = 'Creating mae data',
                 detail = 'This may take a while...', value = 0, {
                     distal.probes <- get.feature.probe(genome = isolate({input$elmerInputGenome}),
                                                        met.platform = isolate({input$elmerInputMetPlatform}))

                     print("Creating MAE")
                     mae <- createMAE(met = met,
                                      exp = exp,
                                      TCGA = isolate({input$elmerInputTCGA}),
                                      filter.probes = distal.probes,
                                      save = FALSE,
                                      met.platform = isolate({input$elmerInputMetPlatform}),
                                      genome = isolate({input$elmerInputGenome}))
                     incProgress(1/2, detail = paste0('MAE created. Saving'))
                     print("Saving MAE")
                     # Define where to save the mae object
                     getPath <- parseDirPath(get.volumes(isolate({input$workingDir})), input$workingDir)
                     if (length(getPath) == 0) getPath <- paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI")
                     filename <- file.path(getPath,isolate({input$maesavefilename}))
                     save(mae,file = filename)
                     incProgress(1/2, detail = paste0('Saving is done'))
                     createAlert(session, "elmerinputmessage", "elmerAlert", title = "MAE created", style =  "success",
                                 content =   paste0("MAE file created: ", filename), append = TRUE)
                 })
    return(mae)

})

observeEvent(input$elmercreatemae, {
    mae <- createMae()
    output$elmerMaeSampleMapping <- DT::renderDataTable({
        return(createTable(as.data.frame(sampleMap(mae))))
    })

    output$elmerMaeSampleMetada <- DT::renderDataTable({
        return(createTable(as.data.frame(colData(mae))))
    })

})
# Input data
elmer.results.data <-  reactive({
    inFile <- input$elmerresultsfile
    if (is.null(inFile) || inFile == 0) return(NULL)
    file  <- as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), inFile)$datapath)

    withProgress(message = 'Loading data',
                 detail = 'This may take a while...', value = 0, {
                     mae.results <- mget(load(file))
                     incProgress(1, detail = "Completed")
                 })
    return(mae.results)
})

maedata <-  reactive({
    inFile <- input$elmermaefile
    if (is.null(inFile) || inFile == 0) return(NULL)
    file  <- as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), inFile)$datapath)

    withProgress(message = 'Loading data',
                 detail = 'This may take a while...', value = 0, {
                     mae <- get(load(file))
                     incProgress(1, detail = "Completed")
                 })
    return(mae)
})

# Updates based on uploaded data
observe({
    results <- elmer.results.data()
    updateSelectizeInput(session, 'scatter.plot.probes', choices = {
        if(!is.null(results)) sort(as.character(results$Sig.probes$probe))
    }, server = TRUE)
    updateSelectizeInput(session, 'scatter.plot.tf', choices = {
        if(!is.null(results)) sort(as.character(rownames(results$TF.meth.cor)))
    }, server = TRUE)
    updateSelectizeInput(session, 'scatter.plot.motif', choices = {
        if(!is.null(results)) sort(as.character(names(results$enriched.motif)))
    }, server = TRUE)
    updateSelectizeInput(session, 'ranking.plot.tf', choices = {
        if(!is.null(results)) sort(as.character(rownames(results$TF.meth.cor)))
    }, server = TRUE)
    updateSelectizeInput(session, 'ranking.plot.motif', choices = {
        if(!is.null(results)) sort(as.character(colnames(results$TF.meth.cor)))
    }, server = TRUE)

})

observeEvent(input$scatter.plot.probes, {
    results <- elmer.results.data()
    updateSelectizeInput(session, 'scatter.plot.genes', choices = {
        if(!is.null(results)) as.character(results$nearGenes[[input$scatter.plot.probes]]$GeneID)
    }, server = TRUE)
})

observe({
    updateSelectizeInput(session, 'schematic.plot.probes', choices = {
        results <- elmer.results.data()
        pair <- results$pair
        if(!is.null(results)){
            as.character(pair$Probe)
        }
    }, server = TRUE)

    updateSelectizeInput(session, 'elmer.heatmap.annotations', choices = {
        results <- elmer.results.data()
        if(!is.null(results)){
            as.character(colnames(colData(results$mae)))
        }
    }, server = TRUE)

})
observeEvent(input$elmermode, {
    updateNumericInput(session, 'elmermetpercentage', value = ifelse(input$elmermode == "Supervised",1,0.2))
    updateNumericInput(session, 'elmergetpairpercentage', value = ifelse(input$elmermode == "Supervised",1,0.4))
    updateNumericInput(session, 'elmergetTFpercentage', value = ifelse(input$elmermode == "Supervised",1,0.4))
})
observe({
    updateSelectizeInput(session, 'schematic.plot.genes', choices = {
        results <- elmer.results.data()
        pair <- results$pair
        if(!is.null(results)){
            genes <-  as.character(pair$GeneID)
            names(genes) <- as.character(paste0(pair$GeneID,"(", pair$Symbol,")"))
            genes
        }
    }, server = TRUE)
})
observeEvent(input$elmerAnalysisBt, {
    closeAlert(session, "elmerAlert")
    getPath <- parseDirPath(get.volumes(isolate({input$workingDir})), input$workingDir)
    if (length(getPath) == 0) getPath <- paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI")
    mae <- maedata()
    if(is.null(mae)){
        closeAlert(session, "elmerAlert")
        createAlert(session, "elmermessage", "elmerAlert", title = "MAE object missing", style =  "danger",
                    content =   "Please upload the MAE object for this plot", append = TRUE)
        return(NULL)
    }



    if(isolate({input$elmergetpairpermu}) > nrow(getMet(mae))){
        closeAlert(session, "elmerAlert")
        createAlert(session, "elmermessage", "elmerAlert", title = "Permutation number problem", style =  "danger",
                    content =   "The number of permutation is higher than the number of probes avaiable, please, reduce it", append = TRUE)
        return(NULL)
    }
    direction <- isolate({input$elmerdirection})
    if(length(direction) == 0){
        closeAlert(session, "elmerAlert")
        createAlert(session, "elmermessage", "elmerAlert", title = "Direction not set", style =  "danger",
                    content =   "Please select at least one direction for the step 1", append = TRUE)
        return(NULL)
    }
    done <- c()
    group1 <- isolate({input$elmermaesubtype})
    group2 <- isolate({input$elmermaesubtype2})
    group.col <- isolate({input$elmermaetype})
    mae <- mae[,colData(mae)[[group.col]] %in% c(group1, group2) ]
    for (j in direction){
        withProgress(message = 'ELMER analysis',
                     detail = paste0('Direction: ',j), value = 0, {
                         print(j)
                         dir.out <- paste0(getPath,"/",isolate({input$elmerresultssavefolder}),"/",j)
                         dir.create(dir.out, recursive = TRUE,showWarnings = FALSE)
                         #--------------------------------------
                         # STEP 3: Analysis                     |
                         #--------------------------------------
                         # Step 3.1: Get diff methylated probes |
                         #--------------------------------------
                         setProgress(value = which(j == direction) * 0.0, message = paste("Step 1-", j, "direction"),
                                     detail = paste("Identify distal probes", j,"methyladed in", group1, "compared to",group2),
                                     session = getDefaultReactiveDomain())

                         diff.dir <- j
                         sig.diff <- isolate({input$elmermetdiff})
                         met.pvalue <- isolate({input$elmermetpvalue})
                         Sig.probes <- get.diff.meth(data = mae,
                                                     sig.dif = sig.diff,
                                                     group.col = group.col,
                                                     group1 = group1,
                                                     group2 = group2,
                                                     cores = isolate({input$elmercores}),
                                                     dir.out = dir.out,
                                                     diff.dir = j,
                                                     mode =  isolate({input$elmermode}),
                                                     pvalue = met.pvalue)
                         if(all(is.na(Sig.probes$probe))){
                             createAlert(session, "elmermessage", "elmerAlert", title = "Error", style =  "error",
                                         content = paste0("No signigicant probes found for ", j ," direction"), append = TRUE)
                             return(NULL)
                         }
                         all.probes <- readr::read_csv(dir(path = dir.out,pattern = "probes.csv",full.names = T,recursive = T,ignore.case = T))
                         #-------------------------------------------------------------
                         # Step 3.2: Identify significant probe-gene pairs            |
                         #-------------------------------------------------------------
                         # Collect nearby 20 genes for Sig.probes
                         setProgress(value = which(j == direction) * 0.1, message = paste("Step 2-", j, "direction"),
                                     detail = paste0("Get Near Genes for ", length(na.omit(Sig.probes$probe))," probes"),
                                     session = getDefaultReactiveDomain())
                         nearGenes <- GetNearGenes(data = mae,
                                                   probes = Sig.probes$probe,
                                                   cores = isolate({input$elmercores}),
                                                   numFlankingGenes = isolate({input$elmergetpairNumGenes}))

                         setProgress(value = which(j == direction) * 0.2, message = paste("Step 3-", j, "direction"),
                                     detail = paste0("Identify putative target genes for differentially methylated distal enhancer probes ", length(na.omit(Sig.probes$probe))," probes"),
                                     session = getDefaultReactiveDomain())

                         pair  <- tryCatch({
                             get.pair(data = mae,
                                      mode =  isolate({input$elmermode}),
                                      nearGenes = nearGenes,
                                      permu.dir = paste0(dir.out,"/permu"),
                                      dir.out = dir.out,
                                      cores = isolate({input$elmercores}),
                                      label = j,
                                      diff.dir = j,
                                      group.col = group.col,
                                      group1 = group1,
                                      group2 = group2,
                                      Pe = isolate({input$elmergetpairpevalue}),
                                      save = TRUE,
                                      raw.pvalue = isolate({input$elmergetpairpvalue}),
                                      permu.size = isolate({input$elmergetpairpermu}),
                                      minSubgroupFrac = isolate({input$elmergetpairpercentage}),
                                      filter.portion = isolate({input$elmergetpairportion}))
                         }, error = function(e) {
                             createAlert(session, "elmermessage", "elmerAlert", title = "Error", style =  "danger",
                                         content = paste0("Error in get.pair function",e), append = TRUE)
                             return(NULL)
                         })


                         if(file.exists(paste0(dir.out,"/getPair.",j,".pairs.significant.csv"))) {
                             Sig.probes.paired <- read.csv(paste0(dir.out,"/getPair.",j,".pairs.significant.csv"),
                                                           stringsAsFactors=FALSE)[,1]
                         } else {
                             createAlert(session, "elmermessage", "elmerAlert", title = "Error", style =  "danger",
                                         content = paste0("No significant pairs were found"), append = TRUE)
                             return(NULL)
                         }
                         #-------------------------------------------------------------
                         # Step 3.3: Motif enrichment analysis on the selected probes |
                         #-------------------------------------------------------------
                         if(length(Sig.probes.paired) > 0 ){
                             #-------------------------------------------------------------
                             # Step 3.3: Motif enrichment analysis on the selected probes |
                             #-------------------------------------------------------------
                             setProgress(value = which(j == direction) * 0.3, message = paste("Step 4-", j, "direction"),
                                         detail = paste0("Identify enriched motifs for the distal enhancer probes which are significantly differentially methylated and linked to putative target gene"),
                                         session = getDefaultReactiveDomain())
                             distal.probe <- get.feature.probe(genome = metadata(mae)$genome,met.platform = metadata(mae)$met.platform)
                             enriched.motif <- get.enriched.motif(data=mae,
                                                                  probes = Sig.probes.paired,
                                                                  dir.out = dir.out,
                                                                  label = j,
                                                                  pvalue = isolate({input$elmergetenrichedmotifFDR}),
                                                                  min.motif.quality = "DS",
                                                                  background.probes = names(distal.probe),
                                                                  lower.OR =  isolate({input$elmergetenrichedmotifLoweOR}),
                                                                  min.incidence = isolate({input$elmergetenrichedmotifMinIncidence}))
                             motif.enrichment <- read.csv(paste0(dir.out,"/getMotif.",j,".motif.enrichment.csv"),
                                                          stringsAsFactors=FALSE)
                             if(length(enriched.motif) > 0){
                                 #-------------------------------------------------------------
                                 # Step 3.4: Identifying regulatory TFs                        |
                                 #-------------------------------------------------------------
                                 print("get.TFs")
                                 setProgress(value = which(j == direction) * 0.4, message = paste("Step 5-", j, "direction"),
                                             detail = paste0(": Identify regulatory TFs whose expression associate with DNA methylation at motifs."), session = getDefaultReactiveDomain())

                                 TF <- get.TFs(data = mae,
                                               mode =  isolate({input$elmermode}),
                                               enriched.motif = enriched.motif,
                                               dir.out = dir.out,
                                               group.col = group.col,
                                               group1 = group1,
                                               group2 = group2,
                                               diff.dir = j,
                                               cores = isolate({input$elmercores}),
                                               label = j,
                                               minSubgroupFrac = isolate({input$elmergetTFpercentage}))
                                 TF.meth.cor <- get(load(paste0(dir.out,"/getTF.",j,".TFs.with.motif.pvalue.rda")))
                                 setProgress(value = which(j == direction) * 0.4, message = paste("Saving results - ", j, "direction"),
                                             detail = paste0("Saving as :", dir.out,"/ELMER_results_",j,".rda"),
                                             session = getDefaultReactiveDomain())

                                 save(mae,
                                      TF,
                                      enriched.motif,
                                      Sig.probes.paired,
                                      pair,
                                      nearGenes,
                                      sig.diff,
                                      met.pvalue,
                                      all.probes,
                                      diff.dir,
                                      Sig.probes,
                                      motif.enrichment,
                                      TF.meth.cor,
                                      group.col,
                                      group1,
                                      group2,
                                      file = paste0(dir.out,"/ELMER_results_",j,".rda"),
                                      compress = "xz")
                                 done <- c(done,j)
                             } else {
                                 createAlert(session, "elmermessage", "elmerAlert", title = "Error", style =  "danger",
                                             content = paste0("No enriched motif was found"), append = TRUE)
                                 return(NULL)
                             }

                         } else {
                             createAlert(session, "elmermessage", "elmerAlert", title = "Error", style =  "danger",
                                         content = paste0("No significant pair gene/probe found"), append = TRUE)
                             return(NULL)
                         }
                         setProgress(value = which(j == direction) * 0.5, message = "Step 5", detail = paste0('Analysis in direction completed: ',j), session = getDefaultReactiveDomain())
                     })
    }
    createAlert(session, "elmermessage", "elmerAlert", title = "ELMER analysis completed", style =  "success",
                content = paste0("ELMER analysis were completed. Results saved in\n",paste0(dir.out,"/ELMER_results_",done,".rda\n", collapse = "")), append = TRUE)
})

elmer.plot <- reactive({
    input$elmerPlotBt
    plot.type <- isolate({input$elmerPlotType})
    mae <- NULL
    results <- elmer.results.data()
    if(!is.null(results)) mae <- results$mae
    closeAlert(session, "elmerAlertResults")

    # Three types:
    # 1 - TF expression vs average DNA methylation
    # 2 - Generate a scatter plot for one probe-gene pair
    # 3 - Generate scatter plots for one probes’ nearby 20 gene expression
    #     vs DNA methylation at this probe
    if(plot.type == "scatter.plot"){
        if(is.null(mae)){
            createAlert(session, "elmermessageresults", "elmerAlertResults", title = "MAE object missing", style =  "danger",
                        content = "Please upload the mae object for this plot", append = TRUE)
            return(NULL)
        }
        # case 1
        plot.by <- isolate({input$scatter.plot.type})
        if(plot.by == "tf"){
            if(is.null(isolate({input$scatter.plot.tf}))){
                closeAlert(session, "elmermessageresults")
                createAlert(session, "elmermessage", "elmerAlertResults", title = "TFs missing", style =  "danger",
                            content = "Please select two TF", append = TRUE)
                return(NULL)
            }

            if(nchar(isolate({input$scatter.plot.tf})) == 0 | length(isolate({input$scatter.plot.tf})) < 2){
                closeAlert(session, "elmerAlertResults")
                createAlert(session, "elmermessageresults", "elmerAlertResults", title = "TFs missing", style =  "danger",
                            content =   "Please select two TF", append = TRUE)
                return(NULL)
            }
            if(nchar(isolate({input$scatter.plot.motif})) == 0){
                closeAlert(session, "elmerAlertResults")
                createAlert(session, "elmermessageresults", "elmerAlertResults", title = "Motif missing", style =  "danger",
                            content =   "Please select a motif", append = TRUE)
                return(NULL)
            }
            p <- scatter.plot(mae,byTF=list(TF=isolate({input$scatter.plot.tf}),
                                            probe=results$enriched.motif[[isolate({input$scatter.plot.motif})]]),
                              category = results$group.col,
                              save = FALSE,
                              lm_line = TRUE)
        } else if(plot.by == "pair") {
            if(nchar(isolate({input$scatter.plot.probes})) == 0){
                closeAlert(session, "elmerAlertResults")
                createAlert(session, "elmermessageresults", "elmerAlertResults", title = "Probe missing", style =  "danger",
                            content =   "Please select a probe", append = TRUE)
                return(NULL)
            }
            if(nchar(isolate({input$scatter.plot.genes})) == 0){
                closeAlert(session, "elmerAlertResults")
                createAlert(session, "elmermessageresults", "elmerAlertResults", title = "Gene missing", style =  "danger",
                            content =   "Please select a gene", append = TRUE)
                return(NULL)
            }

            # case 2
            p <- scatter.plot(mae,byPair=list(probe=isolate({input$scatter.plot.probes}),
                                              gene=c(isolate({input$scatter.plot.genes}))),
                              category=results$group.col, save=FALSE,lm_line=TRUE)
        } else {
            # case 3
            if(nchar(isolate({input$scatter.plot.probes})) == 0){
                createAlert(session, "elmermessageresults", "elmerAlertResults", title = "Probe missing", style =  "danger",
                            content =   "Please select a probe", append = TRUE)
                return(NULL)
            }
            p <- scatter.plot(mae,byProbe=list(probe=isolate({input$scatter.plot.probes}),
                                               numFlankingGenes=isolate({input$scatter.plot.nb.genes})),
                              category = results$group.col,
                              dir.out ="./ELMER.example/Result/LUSC", save=FALSE)
        }
    } else if (plot.type == "schematic.plot") {
        if(is.null(mae)){
            createAlert(session, "elmermessageresults", "elmerAlertResults", title = "mae object missing", style =  "danger",
                        content =   "Please upload the mae object for this plot", append = TRUE)
            return(NULL)
        }
        # Two cases
        # 1 - By probe
        if(isolate({input$schematic.plot.type}) == "probes"){
            if(nchar(isolate({input$schematic.plot.probes})) == 0){
                createAlert(session, "elmermessageresults", "elmerAlertResults", title = "Probe missing", style =  "danger",
                            content =   "Please select a probe", append = TRUE)
                return(NULL)
            }

            p <- schematic.plot(data = mae, pair=results$pair, byProbe=isolate({input$schematic.plot.probes}),save=FALSE)
        } else if(isolate({input$schematic.plot.type}) == "genes"){
            if(nchar(isolate({input$schematic.plot.genes})) == 0){
                createAlert(session, "elmermessageresults", "elmerAlertResults", title = "Gene missing", style =  "success",
                            content =   "Please select a gene", append = TRUE)
                return(NULL)
            }

            # 2 - By genes
            p <- schematic.plot(data = mae, pair=results$pair, byGene=isolate({input$schematic.plot.genes}),save=FALSE)
        }
        p <- recordPlot()
    } else if(plot.type == "motif.enrichment.plot") {
        p <- motif.enrichment.plot(motif.enrichment=results$motif.enrichment,
                                   summary = isolate({input$motif.enrichment.plot.summary}),
                                   significant=list(OR = isolate({input$motif.enrichment.plot.or}),
                                                    lowerOR = isolate({input$motif.enrichment.plot.loweror})),
                                   save=FALSE)
    } else if(plot.type == "ranking.plot"){
        if(nchar(isolate({input$ranking.plot.motif})) == 0){
            createAlert(session, "elmermessageresults", "elmerAlertResults", title = "Motif missing", style =  "success",
                        content =   "Please select a motif", append = TRUE)
            return(NULL)
        }
        label <- isolate({input$ranking.plot.tf})
        if(is.null(label)) {
            label <- createMotifRelevantTfs("family")[isolate({input$ranking.plot.motif})]
        } else {
            label <- list(c(label, unlist(createMotifRelevantTfs("family")[isolate({input$ranking.plot.motif})])))
            names(label) <- isolate({input$ranking.plot.motif})
        }
        gg <- TF.rank.plot(motif.pvalue = results$TF.meth.cor,
                           motif = isolate({input$ranking.plot.motif}),
                           TF.label = label,
                           save = FALSE)
        # names were not fitting in the plot. Reducing the size
        pushViewport(viewport(height = 0.8, width = 0.8))
        p <- grid.draw(gg[[1]])
    } else if(plot.type == "heatmap.plot"){
        if(is.null(mae)){
            createAlert(session, "elmermessageresults", "elmerAlertResults", title = "mae object missing", style =  "danger",
                        content =   "Please upload the mae object for this plot", append = TRUE)
            return(NULL)
        }
        p <- heatmapPairs(data = mae,
                          group.col = results$group.col,
                          group1 = results$group1,
                          annotation.col = isolate({input$elmer.heatmap.annotations}),
                          group2 = results$group2,
                          pairs = results$pair,
                          filename =  NULL)
    }  else if(plot.type == "volcano.plot"){
        if(is.null(mae)){
            createAlert(session, "elmermessageresults", "elmerAlertResults", title = "mae object missing", style =  "danger",
                        content =   "Please upload the mae object for this plot", append = TRUE)
            return(NULL)
        }
        ylab <- ifelse(is.na(diff.dir),
                       " (FDR corrected P-values) [two tailed test]",
                       " (FDR corrected P-values) [one tailed test]")
        p <- TCGAbiolinks:::TCGAVisualize_volcano(x = as.data.frame(all.probes)[,grep("Minus",colnames(all.probes),value = T)],
                                                  y = all.probes$adjust.p,
                                                  title =  paste0("Volcano plot - Probes ",
                                                                  ifelse(is.na(diff.dir),"differently ",diff.dir),
                                                                  "methylated in ", group1, " vs ", group2,"\n"),
                                                  filename = NULL,
                                                  label =  c("Not Significant",
                                                             paste0("Hypermethylated in ",group1),
                                                             paste0("Hypomethylated in ",group1)),
                                                  ylab =  bquote(-Log[10] ~ .(ylab)),
                                                  xlab =  expression(paste(
                                                      "DNA Methylation difference (",beta,"-values)")
                                                  ),
                                                  x.cut = sig.diff,
                                                  y.cut = met.pvalue)
    }
    p
})

observeEvent(input$elmerPlotBt, {
    output$elmer.plot <- renderPlot({
        plot <- elmer.plot()
        if(class(plot) == "list") {
            plot$plot
        } else if(all(class(plot) %in% c("gtable", "gTree", "grob",   "gDesc"))) {
            gridExtra::grid.arrange(plot)
        } else {
            plot
        }
    })
})
observeEvent(input$elmerPlotBt , {
    updateCollapse(session, "collapelmerresults", open = "Plots")
    output$elmerPlot <- renderUI({
        plotOutput("elmer.plot", width = paste0(isolate({input$elmerwidth}), "%"), height = isolate({input$elmerheight}))
    })})

# Table
observeEvent(input$elmerTableType , {
    updateCollapse(session, "collapelmerresults", open = "Results table")
    output$elmerResult <- DT::renderDataTable({
        results <- elmer.results.data()
        if(!is.null(results)){
            if(input$elmerTableType == "tf"){
                createTable(as.data.frame(results$TF))
            } else if(input$elmerTableType == "sigprobes"){
                createTable(as.data.frame(results$Sig.probes))
            } else if(input$elmerTableType == "motif"){
                createTable(as.data.frame(results$motif.enrichment))
            } else if(input$elmerTableType == "pair"){
                createTable(as.data.frame(results$pair))
            }
        }
    })
})


# ELMER INPUT
output$elmerInputSampleMapDownload <- downloadHandler(
    filename = "elmer_example_sample_mapping.tsv",
    content = function(con) {
        met <- elmer.met()
        if(!is.null(met) ) {
            colData <- data.frame(primary = paste0("sample",1:ncol(met)), group = rep("To be filled",ncol(met)))
            write_tsv(colData,con)
        }
    }
)

output$elmerInputcolDataDownload <- downloadHandler(
    filename = "elmer_example_sample_metadata.tsv",
    content = function(con) {
        exp <- elmer.exp()
        met <- elmer.met()
        if(!is.null(exp) & !is.null(met) ) {
            assay <- c(rep("DNA methylation", ncol(met)),
                       rep("Gene expression", ncol(exp)))
            primary <- rep("SampleX", ncol(met) + ncol(exp))
            colname <- c(colnames(met),colnames(exp))
            sampleMap <- data.frame(assay,primary,colname)
            write_tsv(sampleMap,con)
        }
    }
)



output$saveelmerpicture <- downloadHandler(
    filename = function(){input$elmerPlot.filename},
    content = function(file) {
        p <- elmer.plot()
        if(any(class(p) %in% c("gg","ggplot","list"))) {
            if(tools::file_ext(input$elmerPlot.filename) == "png") {
                device <- function(..., width, height) {
                    grDevices::png(..., width = 10, height = 10,
                                   res = 300, units = "in")
                }
            } else if(tools::file_ext(input$elmerPlot.filename) == "pdf") {
                device <- function(..., width, height) {
                    grDevices::pdf(..., width = 10, height = 10)
                }
            } else if(tools::file_ext(input$elmerPlot.filename) == "svg") {
                device <- function(..., width, height) {
                    grDevices::svg(..., width = 10, height = 10)
                }
            }
            if(class(p) == "list") {
                ggsave(file, plot = p$plot, device = device)
            } else {
                ggsave(file, plot = p, device = device)
            }
        } else {
            if(tools::file_ext(file) == "png") {
                grDevices::png(file, width = 10, height = 10,
                               res = 300, units = "in")
            } else if(tools::file_ext(file) == "pdf") {
                grDevices::pdf(file, width = 10, height = 10)
            } else if(tools::file_ext(file) == "svg") {
                grDevices::svg(file, width = 10, height = 10)
            }
            if(class(p) == "recordedplot") {
                print(p)
            } else if(all(class(p) %in% c("gtable", "gTree", "grob",   "gDesc"))) {
                gridExtra::grid.arrange(p)
            } else {
                ComplexHeatmap::draw(p)
            }
            dev.off()
        }
    })
