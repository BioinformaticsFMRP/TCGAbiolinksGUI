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


# mee create
observe({
    data.met <- mee.met()
    updateSelectizeInput(session, 'elmermeetype', choices = {
        if(!is.null(data.met)) as.character(colnames(colData(data.met)))
    }, server = TRUE)
})
observeEvent(input$elmermeetype, {
    data.met <- mee.met()
    updateSelectizeInput(session, 'elmermeesubtype', choices = {
        if(!is.null(data.met)) as.character(unique(colData(data.met)[,input$elmermeetype]))
    }, server = TRUE)
})
observeEvent(input$elmermeetype, {
    data.met <- mee.met()
    updateSelectizeInput(session, 'elmermeesubtype2', choices = {
        if(!is.null(data.met)) as.character(unique(colData(data.met)[,input$elmermeetype]))
    }, server = TRUE)
})

mee.exp <-  reactive({
    inFile <- input$elmerexpfile
    if (is.null(inFile)) return(NULL)
    file  <- as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), inFile)$datapath)

    withProgress(message = 'Loading data',
                 detail = 'This may take a while...', value = 0, {
                     exp <- get(load(file))
                     incProgress(1, detail = "Completed")
                 })
    return(exp)
})
mee.met <-  reactive({
    inFile <- input$elmermetfile
    if (is.null(inFile)) return(NULL)
    file  <- as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), inFile)$datapath)

    withProgress(message = 'Loading data',
                 detail = 'This may take a while...', value = 0, {
                     met <- get(load(file))
                     incProgress(1, detail = "Completed")
                 })

    if(class(met)!= class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
        createAlert(session, "elmermessage", "elmerAlert", title = "Data input error", style =  "danger",
                    content = paste0("Sorry, but I'm expecting a Summarized Experiment object, but I got a: ",
                                     class(met)), append = FALSE)
        return(NULL)
    }

    return(met)
})


observeEvent(input$elmerpreparemee, {
    exp.elmer <- mee.exp()
    met.elmer <- mee.met()
    column <- isolate({input$elmermeetype})
    subtype1 <-  isolate({input$elmermeesubtype})
    subtype2 <-  isolate({input$elmermeesubtype2})

    if(is.null(column)){
        closeAlert(session, "elmerAlert")
        createAlert(session, "elmermessage", "elmerAlert", title = "Type missing", style =  "danger",
                    content =   "Please select the column with type", append = TRUE)
        return(NULL)
    }
    if(is.null(subtype1) | is.null(subtype2)){
        closeAlert(session, "elmerAlert")
        createAlert(session, "elmermessage", "elmerAlert", title = "Subtype missing", style =  "danger",
                    content =   "Please select the two subtypes", append = TRUE)
        return(NULL)
    }
    if(is.null(exp.elmer) | is.null(met.elmer)){
        closeAlert(session, "elmerAlert")
        createAlert(session, "elmermessage", "elmerAlert", title = "Subtype missing", style =  "danger",
                    content =   "Please upload the two summarized Experiment objects", append = TRUE)
        return(NULL)
    }


    # Data: get only samples of subtype1 and subtype2
    exp.elmer <- subset(exp.elmer, select=(colData(exp.elmer)[,column] %in% c(subtype1,subtype2)))
    met.elmer <- subset(met.elmer, select=(colData(met.elmer)[,column] %in%  c(subtype1,subtype2)))

    # Get barcodes for subtype1
    sample.info <-  colData(exp.elmer)
    samples <- sample.info[sample.info[,column] == isolate({input$elmermeesubtype2}),]$barcode
    withProgress(message = 'Creating mee data',
                 detail = 'This may take a while...', value = 0, {

                     exp.elmer <- TCGAprepare_elmer(exp.elmer, platform = "IlluminaHiSeq_RNASeqV2",save = FALSE)
                     incProgress(1/5, detail = paste0('Gene expression matrix prepared'))

                     met.elmer <- TCGAprepare_elmer(met.elmer, platform = "HumanMethylation450",met.na.cut = isolate({input$elmermetnacut}), save = FALSE)
                     incProgress(1/5, detail = paste0('DNA Methylation matrix prepared'))

                     geneAnnot <- txs()
                     geneAnnot$GENEID <- paste0("ID",geneAnnot$GENEID)
                     geneInfo <- promoters(geneAnnot,upstream = 0, downstream = 0)
                     probe <- get.feature.probe()

                     # create mee object, use @ to access the matrices inside the object
                     mee <- fetch.mee(meth = met.elmer, exp = exp.elmer, TCGA = TRUE, probeInfo = probe, geneInfo = geneInfo)
                     incProgress(1/5, detail = paste0('Mee is done'))

                     # Relabel samples in the mee object: subtype1 is control
                     mee@sample$TN[mee@sample$ID %in% substr(samples,1,15)] <- "Control"

                     # Define where to save the mee object
                     getPath <- parseDirPath(get.volumes(isolate({input$workingDir})), input$workingDir)
                     if (length(getPath) == 0) getPath <- paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI")
                     filename <- file.path(getPath,isolate({input$meesavefilename}))
                     save(mee,file = filename)
                     incProgress(2/5, detail = paste0('Saving is done'))
                     createAlert(session, "elmermessage", "elmerAlert", title = "Mee created", style =  "success",
                                 content =   paste0("Mee file created: ", filename), append = TRUE)
                 })
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

meedata <-  reactive({
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
        if(!is.null(results)){
            pair.obj <- fetch.pair(pair=results$pair,
                                   probeInfo =getProbeInfo(results$mee),
                                   geneInfo = getGeneInfo(results$mee))
            as.character(pair.obj@pairInfo$Probe)
        }
    }, server = TRUE)
})
observe({
    updateSelectizeInput(session, 'schematic.plot.genes', choices = {
        results <- elmer.results.data()
        if(!is.null(results)){
            pair.obj <- fetch.pair(pair=results$pair,
                                   probeInfo = getProbeInfo(results$mee),
                                   geneInfo = getGeneInfo(results$mee))

            as.character(pair.obj@pairInfo$GeneID)
        }
    }, server = TRUE)
})
observeEvent(input$elmerAnalysisBt, {
    closeAlert(session, "elmerAlert")
    getPath <- parseDirPath(get.volumes(isolate({input$workingDir})), input$workingDir)
    if (length(getPath) == 0) getPath <- paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI")
    mee <- meedata()
    if(is.null(mee)){
        closeAlert(session, "elmerAlert")
        createAlert(session, "elmermessage", "elmerAlert", title = "Mee object missing", style =  "danger",
                    content =   "Please upload the mee object for this plot", append = TRUE)
        return(NULL)
    }
    if(isolate({input$elmerhyperdir}) & isolate({input$elmerhypodir})) {
        direction <- c("hyper","hypo")
    } else if(isolate({input$elmerhyperdir})) {
        direction <- c("hyper")
    } else {
        direction <- c("hypo")
    }
    done <- c()
    for (j in direction){
        withProgress(message = 'ELMER analysis',
                     detail = paste0('Direction: ',j), value = 0, {
                         print(j)
                         dir.out <- paste0(getPath,"/",isolate({input$elmerresultssavefolder}),"/",j)
                         dir.create(dir.out, recursive = TRUE)
                         #--------------------------------------
                         # STEP 3: Analysis                     |
                         #--------------------------------------
                         # Step 3.1: Get diff methylated probes |
                         #--------------------------------------
                         setProgress(value = which(j == direction) * 0.0, message = "Step 1", detail = paste0(j,": get.diff.meth"), session = getDefaultReactiveDomain())
                         Sig.probes <- get.diff.meth(mee,
                                                     cores=isolate({input$elmercores}),
                                                     dir.out =dir.out,
                                                     diff.dir=j,
                                                     pvalue = isolate({input$elmermetpvalue}))
                         if(all(is.na(Sig.probes$probe))){
                             createAlert(session, "elmermessage", "elmerAlert", title = "Error", style =  "success",
                                         content = paste0("No signigicant probes found"), append = TRUE)
                             return(NULL)
                         }
                         #-------------------------------------------------------------
                         # Step 3.2: Identify significant probe-gene pairs            |
                         #-------------------------------------------------------------
                         # Collect nearby 20 genes for Sig.probes
                         setProgress(value = which(j == direction) * 0.1, message = "Step 2", detail = paste0(j,": GetNearGenes"), session = getDefaultReactiveDomain())
                         nearGenes <- GetNearGenes(TRange=getProbeInfo(mee, probe=Sig.probes$probe),
                                                   cores=isolate({input$elmercores}),
                                                   geneAnnot=getGeneInfo(mee),
                                                   geneNum = isolate({input$elmergetpairNumGenes}))


                         pair  <- tryCatch({
                             get.pair(mee=mee,
                                      probes=Sig.probes$probe,
                                      nearGenes=nearGenes,
                                      permu.dir=paste0(dir.out,"/permu"),
                                      dir.out=dir.out,
                                      cores=isolate({input$elmercores}),
                                      label= j,
                                      permu.size=isolate({input$elmergetpairpermu}),
                                      Pe = isolate({input$elmergetpairpvalue}),
                                      percentage =  isolate({input$elmergetpairpercentage}),
                                      portion = isolate({input$elmergetpairportion}),
                                      diffExp = isolate({input$elmergetpairdiffExp}))
                         }, error = function(e) {
                             createAlert(session, "elmermessage", "elmerAlert", title = "Error", style =  "danger",
                                         content = paste0("Error in get.par function"), append = TRUE)
                             return(NULL)
                         })

                         setProgress(value = which(j == direction) * 0.2, message = "Step 3", detail = paste0(j,": fetch.pair"), session = getDefaultReactiveDomain())
                         Sig.probes.paired <- fetch.pair(pair=pair,
                                                         probeInfo = getProbeInfo(mee),
                                                         geneInfo = getGeneInfo(mee))

                         Sig.probes.paired <- read.csv(paste0(dir.out,"/getPair.",j,".pairs.significant.csv"),
                                                       stringsAsFactors=FALSE)[,1]

                         #-------------------------------------------------------------
                         # Step 3.3: Motif enrichment analysis on the selected probes |
                         #-------------------------------------------------------------
                         if(length(Sig.probes.paired) > 0 ){
                             #-------------------------------------------------------------
                             # Step 3.3: Motif enrichment analysis on the selected probes |
                             #-------------------------------------------------------------
                             setProgress(value = which(j == direction) * 0.3, message = "Step 4", detail = paste0(j,": get.enriched.motif"), session = getDefaultReactiveDomain())
                             probe <- get.feature.probe()
                             enriched.motif <- get.enriched.motif(probes = Sig.probes.paired,
                                                                  dir.out = dir.out,
                                                                  label = j,
                                                                  background.probes = probe$name,
                                                                  lower.OR =  isolate({input$elmergetenrichedmotifLoweOR}),
                                                                  min.incidence = isolate({input$elmergetenrichedmotifMinIncidence}))
                             motif.enrichment <- read.csv(paste0(dir.out,"/getMotif.",j,".motif.enrichment.csv"),
                                                          stringsAsFactors=FALSE)
                             if(length(enriched.motif) > 0){
                                 #-------------------------------------------------------------
                                 # Step 3.4: Identifying regulatory TFs                        |
                                 #-------------------------------------------------------------
                                 print("get.TFs")
                                 setProgress(value = which(j == direction) * 0.4, message = "Step 5", detail = paste0(j,": get.TFs"), session = getDefaultReactiveDomain())

                                 TF <- get.TFs(mee = mee,
                                               enriched.motif = enriched.motif,
                                               dir.out = dir.out,
                                               cores = isolate({input$elmercores}),
                                               label = j,
                                               percentage = isolate({input$elmergetTFpercentage}))
                                 TF.meth.cor <- get(load(paste0(dir.out,"/getTF.",j,".TFs.with.motif.pvalue.rda")))
                                 save(mee, TF, enriched.motif, Sig.probes.paired,
                                      pair, nearGenes, Sig.probes, motif.enrichment, TF.meth.cor,
                                      file=paste0(dir.out,"/ELMER_results_",j,".rda"))
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
                content = paste0("ELMER analysis were completed. Results saved in\n",paste0(dir.out,"/ELMER_results_",done,".rda\n")), append = TRUE)
})

observeEvent(input$elmerPlotBt , {
    output$elmer.plot <- renderPlot({
        plot.type <- isolate({input$elmerPlotType})
        mee <- NULL
        results <- elmer.results.data()
        if(!is.null(results)) mee <- results$mee
        closeAlert(session, "elmerAlert")

        # Three types:
        # 1 - TF expression vs average DNA methylation
        # 2 - Generate a scatter plot for one probe-gene pair
        # 3 - Generate scatter plots for one probes’ nearby 20 gene expression
        #     vs DNA methylation at this probe
        if(plot.type == "scatter.plot"){
            if(is.null(mee)){
                createAlert(session, "elmermessage", "elmerAlert", title = "Mee object missing", style =  "danger",
                            content =   "Please upload the mee object for this plot", append = TRUE)
                return(NULL)
            }
            # case 1
            plot.by <- isolate({input$scatter.plot.type})
            if(plot.by == "tf"){
                if(is.null(isolate({input$scatter.plot.tf}))){
                    closeAlert(session, "elmerAlert")
                    createAlert(session, "elmermessage", "elmerAlert", title = "TFs missing", style =  "danger",
                                content = "Please select two TF", append = TRUE)
                    return(NULL)
                }

                if(nchar(isolate({input$scatter.plot.tf})) == 0 | length(isolate({input$scatter.plot.tf})) < 2){
                    closeAlert(session, "elmerAlert")
                    createAlert(session, "elmermessage", "elmerAlert", title = "TFs missing", style =  "danger",
                                content =   "Please select two TF", append = TRUE)
                    return(NULL)
                }
                if(nchar(isolate({input$scatter.plot.motif})) == 0){
                    closeAlert(session, "elmerAlert")
                    createAlert(session, "elmermessage", "elmerAlert", title = "Motif missing", style =  "danger",
                                content =   "Please select a motif", append = TRUE)
                    return(NULL)
                }
                scatter.plot(mee,byTF=list(TF=isolate({input$scatter.plot.tf}),
                                           probe=results$enriched.motif[[isolate({input$scatter.plot.motif})]]), category="TN",
                             save=FALSE,lm_line=TRUE)
            } else if(plot.by == "pair") {
                if(nchar(isolate({input$scatter.plot.probes})) == 0){
                    closeAlert(session, "elmerAlert")
                    createAlert(session, "elmermessage", "elmerAlert", title = "Probe missing", style =  "danger",
                                content =   "Please select a probe", append = TRUE)
                    return(NULL)
                }
                if(nchar(isolate({input$scatter.plot.genes})) == 0){
                    closeAlert(session, "elmerAlert")
                    createAlert(session, "elmermessage", "elmerAlert", title = "Gene missing", style =  "danger",
                                content =   "Please select a gene", append = TRUE)
                    return(NULL)
                }

                # case 2
                scatter.plot(mee,byPair=list(probe=isolate({input$scatter.plot.probes}),gene=c(isolate({input$scatter.plot.genes}))),
                             category="TN", save=FALSE,lm_line=TRUE)
            } else {
                # case 3
                if(nchar(isolate({input$scatter.plot.probes})) == 0){
                    createAlert(session, "elmermessage", "elmerAlert", title = "Probe missing", style =  "danger",
                                content =   "Please select a probe", append = TRUE)
                    return(NULL)
                }
                scatter.plot(mee,byProbe=list(probe=isolate({input$scatter.plot.probes}),geneNum=isolate({input$scatter.plot.nb.genes})),
                             category="TN", dir.out ="./ELMER.example/Result/LUSC", save=FALSE)
            }
        } else if (plot.type == "schematic.plot") {
            if(is.null(mee)){
                createAlert(session, "elmermessage", "elmerAlert", title = "Mee object missing", style =  "danger",
                            content =   "Please upload the mee object for this plot", append = TRUE)
                return(NULL)
            }
            # Two cases
            # 1 - By probe
            pair.obj <- fetch.pair(pair=results$pair,
                                   probeInfo = getProbeInfo(mee),
                                   geneInfo = getGeneInfo(mee))
            if(isolate({input$schematic.plot.type}) == "probes"){
                if(nchar(isolate({input$schematic.plot.probes})) == 0){
                    createAlert(session, "elmermessage", "elmerAlert", title = "Probe missing", style =  "danger",
                                content =   "Please select a probe", append = TRUE)
                    return(NULL)
                }

                schematic.plot(pair=pair.obj, byProbe=isolate({input$schematic.plot.probes}),save=FALSE)
            } else if(isolate({input$schematic.plot.type}) == "genes"){
                if(nchar(isolate({input$schematic.plot.genes})) == 0){
                    createAlert(session, "elmermessage", "elmerAlert", title = "Gene missing", style =  "success",
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
                createAlert(session, "elmermessage", "elmerAlert", title = "Motif missing", style =  "success",
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
