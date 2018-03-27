#----------------------------------------------------------------------
#   Oncoprint plot
#----------------------------------------------------------------------
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Handling visibility
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observeEvent(input$maftoolsPlot, {
    if(isolate({input$maftoolsPlot}) == "plotmafSummary") {
        shinyjs::show("maftoolsTop")
        shinyjs::hide("maftoolstsb")
        shinyjs::hide("maftoolsgene")
        shinyjs::hide("maftoolsrmNonMutated")
        shinyjs::hide("rainfallPlotps")
    } else  if(isolate({input$maftoolsPlot}) == "oncoplot") {
        shinyjs::show("maftoolsTop")
        shinyjs::hide("maftoolstsb")
        shinyjs::hide("maftoolsgene")
        shinyjs::show("maftoolsrmNonMutated")
        shinyjs::hide("rainfallPlotps")
    } else  if(isolate({input$maftoolsPlot}) == "titv") {
        shinyjs::hide("maftoolsTop")
        shinyjs::hide("maftoolstsb")
        shinyjs::hide("maftoolsgene")
        shinyjs::hide("maftoolsrmNonMutated")
        shinyjs::hide("rainfallPlotps")
    } else  if(isolate({input$maftoolsPlot}) == "rainfallPlot") {
        shinyjs::hide("maftoolsTop")
        shinyjs::show("maftoolstsb")
        shinyjs::hide("maftoolsgene")
        shinyjs::hide("maftoolsrmNonMutated")
        shinyjs::show("rainfallPlotps")
    } else  if(isolate({input$maftoolsPlot}) == "geneCloud") {
        shinyjs::hide("maftoolsTop")
        shinyjs::hide("maftoolstsb")
        shinyjs::hide("maftoolsrmNonMutated")
        shinyjs::hide("rainfallPlotps")
        shinyjs::hide("maftoolsgene")
    } else  if(isolate({input$maftoolsPlot}) == "survivalPlot") {
        shinyjs::hide("maftoolsTop")
        shinyjs::hide("maftoolstsb")
        shinyjs::hide("maftoolsrmNonMutated")
        shinyjs::hide("rainfallPlotps")
        shinyjs::show("maftoolsgene")
    }
})


# The used added a file, now we can show the genes

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# UPDATING FIELDS AFTER DATA INPUT
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# After user uploaded he MAF file we will add the choices availble

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# File selection
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observe({
    shinyFileChoose(input, 'maffile2', roots=get.volumes(input$workingDir), session=session,
                    restrictions=system.file(package='base'), filetypes=c('', "rda",'maf',"csv","maf.gz"))
    shinyFileChoose(input, 'clinicalmaf', roots=get.volumes(input$workingDir), session=session,
                    restrictions=system.file(package='base'), filetypes=c('', "rda","csv"))
})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# DATA INPUT
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
mut.read.maf <-  reactive({
    inFile <- parseFilePaths(get.volumes(isolate({input$workingDir})), input$maffile2)
    if (nrow(inFile) == 0) return(NULL)
    file <- as.character(inFile$datapath)
    withProgress(message = 'Reading MAF file',
                 detail = 'This may take a while...', value = 0, {
                     if(grepl("\\.maf", file)){
                         ret <- read_tsv(file,
                                         comment = "#",
                                         col_types = cols(
                                             Entrez_Gene_Id = col_integer(),
                                             Start_Position = col_integer(),
                                             End_Position = col_integer(),
                                             t_depth = col_integer(),
                                             t_ref_count = col_integer(),
                                             t_alt_count = col_integer(),
                                             n_depth = col_integer(),
                                             ALLELE_NUM = col_integer(),
                                             TRANSCRIPT_STRAND = col_integer(),
                                             PICK = col_integer(),
                                             TSL = col_integer(),
                                             HGVS_OFFSET = col_integer(),
                                             MINIMISED = col_integer()))
                     } else   if(grepl("\\.csv", file)){
                         ret <- read_csv(file)
                     } else {
                         ret <- get(load(file))
                     }
                 })
    return(ret)
})
mut.read.clinical <-  reactive({
    inFile <- parseFilePaths(get.volumes(isolate({input$workingDir})), input$clinicalmaf)
    if (nrow(inFile) == 0) {
        return(NULL)
    } else {
        withProgress(message = 'Reading Clinical file',
                     detail = 'This may take a while...', value = 0, {

                         file <- as.character(inFile$datapath)
                         if(grepl("\\.csv", file)){
                             clinical <- read_csv(file)
                         } else {
                             clinical <- get(load(file))
                         }
                         if(!"Tumor_Sample_Barcode" %in% colnames(clinical)){ # handling TCGA data
                             colnames(clinical)[1] <- "Tumor_Sample_Barcode"
                             clinical$Overall_Survival_Status <- 1
                             clinical$Overall_Survival_Status[which(clinical$vital_status != "dead")] <- 0
                             clinical$time <- clinical$days_to_death
                             clinical$time[is.na(clinical$days_to_death)] <- clinical$days_to_last_follow_up[is.na(clinical$days_to_death)]

                         } else {
                             if(!"Overall_Survival_Status" %in% colnames(clinical)){
                                 createAlert(session, "maftoolsmessage", "maftoolsAlert", title = "Error", style =  "danger",
                                             content = "Please we did not find column Overall_Survival_Status in the file (1 = dead, 0 alive)", append = TRUE)
                                 return(NULL)
                             }
                             if(!"time" %in% colnames(clinical)){
                                 createAlert(session, "maftoolsmessage", "maftoolsAlert", title = "Error", style =  "danger",
                                             content = "Please we did not find column time in the file (days to death)", append = TRUE)
                                 return(NULL)
                             }
                         }
                         incProgress(1, detail = paste("Done"))
                     })
    }
    return(clinical)
})
mut <-  reactive({
    maf <- mut.read.maf()
    if(is.null(maf)) return(NULL)
    clinical <- mut.read.clinical()
    if(is.null(clinical)) {
        ret <- read.maf(maf = maf)
    } else {
        ret <- read.maf(maf = maf, clinicalData = clinical,isTCGA = T)
    }
    return(ret)
})

observe({
    data <- mut()
    updateSelectizeInput(session, 'maftoolstsb', choices = {
        if(!is.null(data)) as.character(unique(getSampleSummary(data)$Tumor_Sample_Barcode))
    }, server = TRUE)
    updateSelectizeInput(session, 'maftoolsgene', choices = {
        if(!is.null(data)) as.character(unique(getGeneSummary(data)$Hugo_Symbol))
    }, server = TRUE)

})
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Plot
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
maftools.plot <-  reactive({
    input$maftoolsPlotBt
    mut <- isolate({mut()})

    if(is.null(mut)){
        createAlert(session, "maftoolsmessage", "maftoolsAlert", title = "Error", style =  "danger",
                    content = "Please select a file", append = TRUE)
        return(NULL)
    }
    withProgress(message = 'Creating plot',
                 detail = 'This may take a while...', value = 0, {

                     if(isolate({input$maftoolsPlot}) == "plotmafSummary") {
                         p <- plotmafSummary(maf = mut,
                                             top = isolate({input$maftoolsTop}),
                                             rmOutlier = TRUE,
                                             addStat = 'median',
                                             dashboard = TRUE)
                     } else  if(isolate({input$maftoolsPlot}) == "oncoplot") {
                         p <- oncoplot(maf = mut, top = isolate({input$maftoolsTop}),
                                       removeNonMutated = isolate({input$maftoolsrmNonMutated}))
                     } else  if(isolate({input$maftoolsPlot}) == "titv") {
                         titv <- titv(maf = mut, plot = FALSE, useSyn = TRUE)
                         p <- plotTiTv(res = titv)
                         p <- recordPlot()
                     } else  if(isolate({input$maftoolsPlot}) == "rainfallPlot") {
                         if(str_length(isolate({input$maftoolstsb})) == 0) {
                             tsb <- NULL
                         } else {
                             tsb <- isolate({input$maftoolstsb})
                         }
                         p <- rainfallPlot(maf = mut,
                                           detectChangePoints = TRUE,
                                           tsb = tsb,
                                           fontSize = 12,
                                           pointSize = isolate({input$rainfallPlotps}) )
                     } else  if(isolate({input$maftoolsPlot}) == "geneCloud") {
                         p <- geneCloud(input = mut, minMut = 3)
                         p <- recordPlot()
                     } else  if(isolate({input$maftoolsPlot}) == "survivalPlot") {
                         genes <- isolate({input$maftoolsgene})
                         p <- mafSurvival(maf = mut,
                                          genes = genes,
                                          fn = genes,
                                          groupNames = c(paste(genes, "Mutant"),"WT"),
                                          time = 'time',
                                          Status = 'Overall_Survival_Status',
                                          isTCGA = TRUE)
                     }
                 })
    p
})

observeEvent(input$maftoolsPlotBt, {
    closeAlert(session,"maftoolsAlert")
    output$maftoolsploting <- renderPlot({
        maftools.plot()
    })
})

observeEvent(input$maftoolsPlotBt , {
    updateCollapse(session, "collapsemaftools", open = "Plots")
    output$maftoolsplot <- renderUI({
        plotOutput("maftoolsploting", height = "600px")
    })
})

output$savemafpicture <- downloadHandler(
    filename = function(){input$mafPlot.filename},
    content = function(file) {
        p <- maftools.plot()
        if(any(class(p) %in% c("gg","ggplot","list"))) {
            if(tools::file_ext(input$mafPlot.filename) == "png") {
                device <- function(..., width, height) {
                    grDevices::png(..., width = 10, height = 10,
                                   res = 300, units = "in")
                }
            } else if(tools::file_ext(input$mafPlot.filename) == "pdf") {
                device <- function(..., width, height) {
                    grDevices::pdf(..., width = 10, height = 10)
                }
            } else if(tools::file_ext(input$mafPlot.filename) == "svg") {
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
            } else {
                ComplexHeatmap::draw(p)
            }
            dev.off()
        }
    })

# Table
observeEvent(input$maftoolsTableBt , {
    updateCollapse(session, "collapsemaftools", open = "Survival log rank test table")
    output$maftoolsSurvivalTbl <- DT::renderDataTable({
        mut <- isolate({mut()})
        if(!is.null(mut)){
            if(nrow(getClinicalData(mut)) > 0) {
                genes <- getGeneSummary(mut)$Hugo_Symbol[1:isolate({input$maftoolsTopMutated})]
                tbl <- plyr::adply(genes,1, function(x){
                    surv = mafSurvival(maf = mut,
                                       genes = x,
                                       fn = x,
                                       time = 'time',
                                       Status = 'Overall_Survival_Status',
                                       isTCGA = TRUE)
                    pvalue <- surv$labels$subtitle
                    pvalue
                },.progress = "text")
                tbl$X1 <- genes
                colnames(tbl) <- c("gene","pvalue")
                tbl$pvalue <- gsub("P-value: ","",tbl$pvalue )
                tbl$pvalue <- as.numeric(tbl$pvalue)
                tbl <- tbl[order(tbl$pvalue),]
                return(createTable(tbl))
            }
        }
    })
})


