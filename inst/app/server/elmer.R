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
        shinyjs::hide("ranking.plot.tf")
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
        shinyjs::hide("scatter.plot.genes")
        shinyjs::hide("scatter.plot.probes")
        shinyjs::hide("scatter.plot.nb.genes")
        shinyjs::show("schematic.plot.type")
        shinyjs::hide("ranking.plot.motif")
        shinyjs::hide("ranking.plot.tf")
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
        shinyjs::hide("scatter.plot.probes")
        shinyjs::hide("scatter.plot.nb.genes")
        shinyjs::hide("schematic.plot.type")
        shinyjs::hide("schematic.plot.genes")
        shinyjs::hide("schematic.plot.probes")
        shinyjs::show("ranking.plot.motif")
        shinyjs::show("ranking.plot.tf")
    } else if(type =="motif.enrichment.plot"){
        shinyjs::hide("scatter.plot.type")
        shinyjs::hide("scatter.plot.tf")
        shinyjs::hide("scatter.plot.motif")
        shinyjs::hide("scatter.plot.genes")
        shinyjs::hide("scatter.plot.probes")
        shinyjs::hide("scatter.plot.nb.genes")
        shinyjs::hide("schematic.plot.type")
        shinyjs::hide("schematic.plot.genes")
        shinyjs::hide("schematic.plot.probes")
        shinyjs::hide("ranking.plot.motif")
        shinyjs::hide("ranking.plot.tf")
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
    shinyFileChoose(input, 'elmermeefile', roots=get.volumes(input$workingDir), session=session,
                    restrictions=system.file(package='base'), filetypes=c('', 'rda'))
    shinyFileChoose(input, 'elmerresultsfile', roots=get.volumes(input$workingDir), session=session,
                    restrictions=system.file(package='base'), filetypes=c('', 'rda'))
    shinyFileChoose(input, 'elmermetfile', roots=get.volumes(input$workingDir), session=session,
                    restrictions=system.file(package='base'), filetypes=c('', 'rda'))
    shinyFileChoose(input, 'elmerexpfile', roots=get.volumes(input$workingDir), session=session,
                    restrictions=system.file(package='base'), filetypes=c('', 'rda'))
    updateTextInput(session, "meesavefilename",
                    value = paste0("mee",
                                   paste(input$elmermeetype,
                                         gsub("[[:punct:]]| ","_",input$elmermeesubtype),
                                         gsub("[[:punct:]]| ","_",input$elmermeesubtype2),
                                         sep = "_"),".rda"))
})


# mae create
observe({
    data <- maedata()
    updateSelectizeInput(session, 'elmermeetype', choices = {
        if(!is.null(data)) as.character(colnames(pData(data)))
    }, server = TRUE)
})
observeEvent(input$elmermeetype, {
    data <- maedata()
    updateSelectizeInput(session, 'elmermeesubtype', choices = {
        if(!is.null(data)) as.character(unique(pData(data)[,input$elmermeetype]))
    }, server = TRUE)
})
observeEvent(input$elmermeetype, {
    data <- maedata()
    updateSelectizeInput(session, 'elmermeesubtype2', choices = {
        if(!is.null(data)) as.character(unique(pData(data)[,input$elmermeetype]))
    }, server = TRUE)
})

observe({
    input$elmermeesubtype2
    input$elmermeesubtype
    group1 <- if_else(str_length(isolate({input$elmermeesubtype})) > 0,isolate({input$elmermeesubtype}), "group 1")
    group2 <- if_else(str_length(isolate({input$elmermeesubtype2})) > 0,isolate({input$elmermeesubtype2}), "group 2")
    choices <- c("hyper","hypo")
    names(choices) <-  c(paste("Probes hypermethylated in", group1,"compared to", group2),
                         paste("Probes hypomethylated in", group1,"compared to", group2))
    updateSelectizeInput(session, 'elmerdirection', choices = {
          choices
    }, server = TRUE)
})

elmer.exp <-  reactive({
    inFile <- input$elmerexpfile
    if (is.null(inFile)) return(NULL)
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
    if (is.null(inFile)) return(NULL)
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

    withProgress(message = 'Creating mee data',
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
                     # Define where to save the mee object
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
    output$elmerMaeSampleMapping <- renderDataTable({
        return(as.data.frame(sampleMap(mae)))
    },
    options = list(pageLength = 10,
                   scrollX = TRUE,
                   jQueryUI = TRUE,
                   pagingType = "full",
                   lengthMenu = list(c(10, 20, -1), c('10', '20', 'All')),
                   language.emptyTable = "No results found",
                   "dom" = 'T<"clear">lfrtip',
                   "oTableTools" = list(
                       "sSelectedClass" = "selected",
                       "sRowSelect" = "os",
                       "sSwfPath" = paste0("//cdn.datatables.net/tabletools/2.2.4/swf/copy_csv_xls.swf"),
                       "aButtons" = list(
                           list("sExtends" = "collection",
                                "sButtonText" = "Save",
                                "aButtons" = c("csv","xls")
                           )
                       )
                   )
    ), callback = "function(table) {table.on('click.dt', 'tr', function() {Shiny.onInputChange('allRows',table.rows('.selected').data().toArray());});}"
    )

    output$elmerMaeSampleMetada <- renderDataTable({
        return(as.data.frame(pData(mae)))
    },
    options = list(pageLength = 10,
                   scrollX = TRUE,
                   jQueryUI = TRUE,
                   pagingType = "full",
                   lengthMenu = list(c(10, 20, -1), c('10', '20', 'All')),
                   language.emptyTable = "No results found",
                   "dom" = 'T<"clear">lfrtip',
                   "oTableTools" = list(
                       "sSelectedClass" = "selected",
                       "sRowSelect" = "os",
                       "sSwfPath" = paste0("//cdn.datatables.net/tabletools/2.2.4/swf/copy_csv_xls.swf"),
                       "aButtons" = list(
                           list("sExtends" = "collection",
                                "sButtonText" = "Save",
                                "aButtons" = c("csv","xls")
                           )
                       )
                   )
    ), callback = "function(table) {table.on('click.dt', 'tr', function() {Shiny.onInputChange('allRows',table.rows('.selected').data().toArray());});}"
    )

})
# Input data
elmer.results.data <-  reactive({
    inFile <- input$elmerresultsfile
    if (is.null(inFile)) return(NULL)
    file  <- as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), inFile)$datapath)

    withProgress(message = 'Loading data',
                 detail = 'This may take a while...', value = 0, {
                     mee.results <- mget(load(file))
                     incProgress(1, detail = "Completed")
                 })
    return(mee.results)
})

maedata <-  reactive({
    inFile <- input$elmermeefile
    if (is.null(inFile)) return(NULL)
    file  <- as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), inFile)$datapath)

    withProgress(message = 'Loading data',
                 detail = 'This may take a while...', value = 0, {
                     mee <- get(load(file))
                     incProgress(1, detail = "Completed")
                 })
    return(mee)
})

# Updates based on uploaded data
observe({
    results <- elmer.results.data()
    updateSelectizeInput(session, 'scatter.plot.probes', choices = {
        if(!is.null(results)) as.character(results$Sig.probes$probe)
    }, server = TRUE)
    updateSelectizeInput(session, 'scatter.plot.tf', choices = {
        if(!is.null(results)) as.character(rownames(results$TF.meth.cor))
    }, server = TRUE)
    updateSelectizeInput(session, 'scatter.plot.motif', choices = {
        if(!is.null(results)) as.character(names(results$enriched.motif))
    }, server = TRUE)
    updateSelectizeInput(session, 'ranking.plot.tf', choices = {
        if(!is.null(results)) as.character(rownames(results$TF.meth.cor))
    }, server = TRUE)
    updateSelectizeInput(session, 'ranking.plot.motif', choices = {
        if(!is.null(results)) as.character(colnames(results$TF.meth.cor))
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
    done <- c()
    group1 <- isolate({input$elmermeesubtype})
    group2 <- isolate({input$elmermeesubtype2})
    group.col <- isolate({input$elmermeetype})
    mae <- mae[,pData(mae)[[group.col]] %in% c(group1, group2) ]
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
                         Sig.probes <- get.diff.meth(mae,
                                                     group.col = group.col,
                                                     group1 = group1,
                                                     group2 = group2,
                                                     cores = isolate({input$elmercores}),
                                                     dir.out = dir.out,
                                                     diff.dir = j,
                                                     pvalue = isolate({input$elmermetpvalue}))
                         if(all(is.na(Sig.probes$probe))){
                             createAlert(session, "elmermessage", "elmerAlert", title = "Error", style =  "error",
                                         content = paste0("No signigicant probes found for ", j ," direction"), append = TRUE)
                             return(NULL)
                         }
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
                                      nearGenes = nearGenes,
                                      permu.dir = paste0(dir.out,"/permu"),
                                      dir.out = dir.out,
                                      cores = isolate({input$elmercores}),
                                      label = j,
                                      group.col = group.col,
                                      group1 = group1,
                                      group2 = group2,
                                      Pe = isolate({input$elmergetpairpvalue}),
                                      save = TRUE,
                                      pvalue = isolate({input$elmergetpairpvalue}),
                                      permu.size=isolate({input$elmergetpairpermu}),
                                      minSubgroupFrac = isolate({input$elmergetpairpercentage}),
                                      filter.portion = isolate({input$elmergetpairportion}),
                                      diffExp = isolate({input$elmergetpairdiffExp}))
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
                                               enriched.motif = enriched.motif,
                                               dir.out = dir.out,
                                               group.col = group.col,
                                               group1 = group1,
                                               group2 = group2,
                                               cores = isolate({input$elmercores}),
                                               label = j,
                                               minSubgroupFrac = isolate({input$elmergetTFpercentage}))
                                 TF.meth.cor <- get(load(paste0(dir.out,"/getTF.",j,".TFs.with.motif.pvalue.rda")))
                                 save(mae, TF, enriched.motif, Sig.probes.paired,
                                      pair, nearGenes, Sig.probes, motif.enrichment, TF.meth.cor,
                                      group.col, group1, group2,
                                      file=paste0(dir.out,"/ELMER_results_",j,".rda"),
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

observeEvent(input$elmerPlotBt , {
    output$elmer.plot <- renderPlot({
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
                scatter.plot(mae,byTF=list(TF=isolate({input$scatter.plot.tf}),
                                           probe=results$enriched.motif[[isolate({input$scatter.plot.motif})]]), category="TN",
                             save=FALSE,lm_line=TRUE)
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
                scatter.plot(mee,byPair=list(probe=isolate({input$scatter.plot.probes}),gene=c(isolate({input$scatter.plot.genes}))),
                             category="TN", save=FALSE,lm_line=TRUE)
            } else {
                # case 3
                if(nchar(isolate({input$scatter.plot.probes})) == 0){
                    createAlert(session, "elmermessageresults", "elmerAlertResults", title = "Probe missing", style =  "danger",
                                content =   "Please select a probe", append = TRUE)
                    return(NULL)
                }
                scatter.plot(mee,byProbe=list(probe=isolate({input$scatter.plot.probes}),geneNum=isolate({input$scatter.plot.nb.genes})),
                             category=results$group.col, dir.out ="./ELMER.example/Result/LUSC", save=FALSE)
            }
        } else if (plot.type == "schematic.plot") {
            if(is.null(mee)){
                createAlert(session, "elmermessageresults", "elmerAlertResults", title = "Mee object missing", style =  "danger",
                            content =   "Please upload the mee object for this plot", append = TRUE)
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

                schematic.plot(pair=pair.obj, byProbe=isolate({input$schematic.plot.probes}),save=FALSE)
            } else if(isolate({input$schematic.plot.type}) == "genes"){
                if(nchar(isolate({input$schematic.plot.genes})) == 0){
                    createAlert(session, "elmermessageresults", "elmerAlertResults", title = "Gene missing", style =  "success",
                                content =   "Please select a gene", append = TRUE)
                    return(NULL)
                }

                # 2 - By genes
                schematic.plot(pair=pair.obj, byGene=isolate({input$schematic.plot.genes}),save=FALSE)
            }
        } else if(plot.type == "motif.enrichment.plot") {
            motif.enrichment.plot(motif.enrichment=results$motif.enrichment,
                                  #significant=list(OR=1.3,lowerOR=1.3),
                                  save=FALSE)
        } else if(plot.type == "ranking.plot"){
            if(nchar(isolate({input$ranking.plot.motif})) == 0){
                createAlert(session, "elmermessageresults", "elmerAlertResults", title = "Motif missing", style =  "success",
                            content =   "Please select a motif", append = TRUE)
                return(NULL)
            }
            label <- list(isolate({input$ranking.plot.tf}))
            names(label) <- isolate({input$ranking.plot.motif})
            gg <- TF.rank.plot(motif.pvalue=results$TF.meth.cor,
                               motif=isolate({input$ranking.plot.motif}),
                               TF.label=label,
                               save=FALSE)
            # names were not fitting in the plot. Reducing the size
            pushViewport(viewport(height=0.8,width=0.8))
            grid.draw(gg[[1]])
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
    output$elmerResult <- renderDataTable({
        results <- elmer.results.data()
        if(!is.null(results)){
            if(input$elmerTableType == "tf"){
                as.data.frame(results$TF)
            } else if(input$elmerTableType == "sigprobes"){
                as.data.frame(results$Sig.probes)
            } else if(input$elmerTableType == "motif"){
                as.data.frame(results$motif.enrichment)
            } else if(input$elmerTableType == "pair"){
                as.data.frame(results$pair)
            }
        }
    },
    options = list(pageLength = 10,
                   scrollX = TRUE,
                   jQueryUI = TRUE,
                   pagingType = "full",
                   lengthMenu = list(c(10, 20, -1), c('10', '20', 'All')),
                   language.emptyTable = "No results found",
                   "dom" = 'T<"clear">lfrtip',
                   "oTableTools" = list(
                       "sSelectedClass" = "selected",
                       "sRowSelect" = "os",
                       "sSwfPath" = paste0("//cdn.datatables.net/tabletools/2.2.4/swf/copy_csv_xls.swf"),
                       "aButtons" = list(
                           list("sExtends" = "collection",
                                "sButtonText" = "Save",
                                "aButtons" = c("csv","xls")
                           )
                       )
                   )
    ), callback = "function(table) {table.on('click.dt', 'tr', function() {Shiny.onInputChange('allRows',table.rows('.selected').data().toArray());});}"
    )
})


# ELMER INPUT
output$elmerInputSampleMapDownload <- downloadHandler(
    filename = "elmer_example_sample_mapping.tsv",
    content = function(con) {
        met <- elmer.met()
        if(!is.null(met) ) {
            pData <- data.frame(primary = paste0("sample",1:ncol(met)), group = rep("To be filled",ncol(met)))
            write_tsv(pData,con)
        }
    }
)

output$elmerInputpDataDownload <- downloadHandler(
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

