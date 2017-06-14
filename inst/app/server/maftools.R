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
        shinyjs::hide("maftoolsrmNonMutated")
        shinyjs::hide("rainfallPlotps")
    } else  if(isolate({input$maftoolsPlot}) == "oncoplot") {
        shinyjs::show("maftoolsTop")
        shinyjs::hide("maftoolstsb")
        shinyjs::show("maftoolsrmNonMutated")
        shinyjs::hide("rainfallPlotps")
    } else  if(isolate({input$maftoolsPlot}) == "titv") {
        shinyjs::hide("maftoolsTop")
        shinyjs::hide("maftoolstsb")
        shinyjs::hide("maftoolsrmNonMutated")
        shinyjs::hide("rainfallPlotps")
    } else  if(isolate({input$maftoolsPlot}) == "rainfallPlot") {
        shinyjs::hide("maftoolsTop")
        shinyjs::show("maftoolstsb")
        shinyjs::hide("maftoolsrmNonMutated")
        shinyjs::show("rainfallPlotps")
    } else  if(isolate({input$maftoolsPlot}) == "geneCloud") {
        shinyjs::hide("maftoolsTop")
        shinyjs::hide("maftoolstsb")
        shinyjs::hide("maftoolsrmNonMutated")
        shinyjs::hide("rainfallPlotps")
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
})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# DATA INPUT
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
mut <-  reactive({
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
                         ret <-  get(load(file))
                     }
                     incProgress(1, detail = paste("Done"))
                 })
    return(ret)
})

observe({
    data <- mut()
    updateSelectizeInput(session, 'maftoolstsb', choices = {
        if(!is.null(data)) as.character(unique(data$Tumor_Sample_Barcode))
    }, server = TRUE)
})
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Plot
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observeEvent(input$maftoolsPlotBt , {
    closeAlert(session,"maftoolsAlert")
    output$maftoolsploting <- renderPlot({
        mut <- isolate({mut()})

        if(is.null(mut)){
            createAlert(session, "oncomessage", "oncoAlert", title = "Error", style =  "danger",
                        content = "Please select a file", append = TRUE)
            return(NULL)
        }
        mut = read.maf(maf = mut, removeSilent = TRUE, useAll = FALSE)
        withProgress(message = 'Creating plot',
                     detail = 'This may take a while...', value = 0, {

                         if(isolate({input$maftoolsPlot}) == "plotmafSummary") {
                             plotmafSummary(maf = mut,top = isolate({input$maftoolsTop}),
                                            rmOutlier = TRUE,
                                            addStat = 'median', dashboard = TRUE)
                         } else  if(isolate({input$maftoolsPlot}) == "oncoplot") {
                             oncoplot(maf = mut, top = isolate({input$maftoolsTop}),
                                      removeNonMutated = isolate({input$maftoolsrmNonMutated}))
                         } else  if(isolate({input$maftoolsPlot}) == "titv") {
                             titv <- titv(maf = mut, plot = FALSE, useSyn = TRUE)
                             plotTiTv(res = titv)
                         } else  if(isolate({input$maftoolsPlot}) == "rainfallPlot") {
                             if(str_length(isolate({input$maftoolstsb})) == 0) {
                                 tsb <- NULL
                             } else {
                                 tsb <- isolate({input$maftoolstsb})
                             }
                             rainfallPlot(maf = mut,
                                          detectChangePoints = TRUE,
                                          tsb = tsb,
                                          fontSize = 12,
                                          pointSize = isolate({input$rainfallPlotps}) )
                         } else  if(isolate({input$maftoolsPlot}) == "geneCloud") {
                             geneCloud(input = mut, minMut = 3)
                         }
                     })
    })
})

observeEvent(input$maftoolsPlotBt , {
    updateCollapse(session, "collapsemaftools", open = "Plots")
    output$maftoolsplot <- renderUI({
        plotOutput("maftoolsploting", height = "600px")
    })
})
