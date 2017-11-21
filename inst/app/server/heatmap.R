##----------------------------------------------------------------------
#                             Heatmap
##----------------------------------------------------------------------
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Handling visibility
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observe({
    if(input$heatmapTypeInputRb == "met") {
        shinyjs::show("heatmapProbesInputRb")
        shinyjs::hide("heatmapGenesInputRb")
        shinyjs::hide("heatmapGenesTextArea")
        shinyjs::hide("heatmagenes")
        shinyjs::hide("heatmap.upGenesCb")
        shinyjs::hide("heatmap.downGenewsCb")
        if(input$heatmapProbesInputRb == "text"){
            shinyjs::show("heatmapProbesTextArea")
            shinyjs::hide("heatmap.hypoprobesCb")
            shinyjs::hide("heatmaprobes")
            shinyjs::hide("heatmap.hyperprobesCb")
        } else if(input$heatmapProbesInputRb == "Selection"){
            shinyjs::hide("heatmapProbesTextArea")
            shinyjs::hide("heatmap.hypoprobesCb")
            shinyjs::hide("heatmap.hyperprobesCb")
            shinyjs::show("heatmaprobes")
        } else {
            shinyjs::hide("heatmapProbesTextArea")
            shinyjs::show("heatmap.hypoprobesCb")
            shinyjs::hide("heatmaprobes")
            shinyjs::show("heatmap.hyperprobesCb")
        }
    } else if(input$heatmapTypeInputRb == "exp") {
        shinyjs::show("heatmapGenesInputRb")
        shinyjs::hide("heatmapProbesTextArea")
        shinyjs::hide("heatmap.hypoprobesCb")
        shinyjs::hide("heatmaprobes")
        shinyjs::hide("heatmap.hyperprobesCb")
        shinyjs::hide("heatmapProbesInputRb")
        if(input$heatmapGenesInputRb == "text"){
            shinyjs::show("heatmapGenesTextArea")
            shinyjs::hide("heatmap.upGenesCb")
            shinyjs::hide("heatmap.downGenewsCb")
            shinyjs::hide("heatmagenes")
        } else if(input$heatmapGenesInputRb == "Selection"){
            shinyjs::hide("heatmapGenesTextArea")
            shinyjs::hide("heatmap.upGenesCb")
            shinyjs::hide("heatmap.downGenewsCb")
            shinyjs::show("heatmagenes")
        } else {
            shinyjs::hide("heatmapGenesTextArea")
            shinyjs::show("heatmap.upGenesCb")
            shinyjs::show("heatmap.downGenewsCb")
            shinyjs::hide("heatmagenes")
        }
    }
})
observeEvent(input$heatmap.colorsCb, {
    if(input$heatmap.colorsCb){
        shinyjs::show("heatmapcolMax")
        shinyjs::show("heatmapcolMid")
        shinyjs::show("heatmapcolMin")
    } else {
        shinyjs::hide("heatmapcolMax")
        shinyjs::hide("heatmapcolMid")
        shinyjs::hide("heatmapcolMin")
    }
})
observeEvent(input$heatmap.extremesCb, {
    if(input$heatmap.extremesCb) {
        shinyjs::show("heatmapExtremeMax")
        shinyjs::show("heatmapExtremeMid")
        shinyjs::show("heatmapExtremeMin")
    } else {
        shinyjs::hide("heatmapExtremeMax")
        shinyjs::hide("heatmapExtremeMid")
        shinyjs::hide("heatmapExtremeMin")
    }
})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# File selection
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observe({
    shinyFileChoose(input, 'heatmapfile', roots=get.volumes(input$workingDir), session=session,
                    restrictions=system.file(package='base'),filetypes=c('', 'rda'))
    shinyFileChoose(input, 'heatmapresultsfile', roots=get.volumes(input$workingDir), session=session,
                    restrictions=system.file(package='base'),filetypes=c('', 'csv'))
})
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# UPDATING FIELDS AFTER DATA INPUT
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observe({
    if((input$heatmap.sortCb)){
        updateCheckboxInput(session, "heatmap.clustercol",  value = FALSE)
    }
})

observe({
    updateSelectizeInput(session, 'heatmapSortCol', choices = {
        if (class(heatmapdata()) == class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
            if (!is.null(heatmapdata()) & !is.null(input$colmetadataheatmap))
                as.character(input$colmetadataheatmap)
        }}, server = TRUE)
})

# If the user ulpload a DEA results, the type of heatmap should be changed automatically to avoid
# errors. Should this be visible to the user ?
observe({
    if(!is.null(input$heatmapresultsfile)){
        file  <- basename(as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), input$heatmapresultsfile)$datapath))
        selected <- "met"
        if(grepl("DEA",file))  selected <- "exp"
        updateRadioButtons(session, "heatmapTypeInputRb", selected = selected)
        updateTextInput(session,"heatmapLabel", value = ifelse(selected == "met","DNA methylation level","Gene expression level"))
    }
})

observe({
    data <- heatmapdata()
    updateSelectizeInput(session, 'colmetadataheatmap', choices = {
        if(!is.null(data)) as.character(colnames(colData(data)))
    }, server = TRUE)
})
observe({
    data <- heatmapdata()
    updateSelectizeInput(session, 'rowmetadataheatmap', choices = {
        if(!is.null(data)) as.character(colnames(values(data)))
    }, server = TRUE)
    updateSelectizeInput(session, 'heatmaprobes', choices = as.character(rownames(data)), server = TRUE)
    updateSelectizeInput(session, 'heatmagenes', choices = as.character(rownames(data)), server = TRUE)
})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# DATA INPUT
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
heatmapresultdata <-  reactive({
    inFile <- input$heatmapresultsfile
    if (is.null(inFile)) return(NULL)
    file  <- as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), inFile)$datapath)
    # verify if the file is a csv
    ext <- tools::file_ext(file)
    if(ext != "csv"){
        createAlert(session, "dmrmessage", "dmrAlert", title = "Data input error", style =  "danger",
                    content = paste0("Sorry, but I'm expecting a csv file, but I got a: ",
                                     ext), append = FALSE)
        return(NULL)
    }

    withProgress(message = 'Loading data',
                 detail = 'This may take a while...', value = 0, {
                     df <- as.data.frame(read_csv(file)); rownames(df) <- df$X1; df$X1 <- NULL;
                     incProgress(1, detail = "Completed")
                 })
    return(df)
})

heatmapdata <-  reactive({
    inFile <- input$heatmapfile
    if (is.null(inFile)) return(NULL)
    file  <- as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), input$heatmapfile)$datapath)

    withProgress(message = 'Loading data',
                 detail = 'This may take a while...', value = 0, {
                     result.file <- gsub(".rda","_results.rda",file)
                     if(file.exists(result.file)) {
                         se <- get(load(result.file))
                     } else {
                         se <- get(load(file))
                     }
                     incProgress(1, detail = "Completed")
                 })
    if(class(se)!= class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
        createAlert(session, "heatmapmessage", "heatmapAlert", title = "Data input error", style =  "danger",
                    content = paste0("Sorry, but I'm expecting a Summarized Experiment object, but I got a: ",
                                     class(se)), append = FALSE)
        return(NULL)
    }
    return(se)

})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Plot
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
heatmap.plot <-  reactive({
    input$heatmapPlotBt
    file  <- basename(as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), input$heatmapresultsfile)$datapath))
    if(length(file) > 0){
        file <- unlist(str_split(file,"_"))
        group1 <- file[4]
        group2 <- file[5]
    }
    data <- isolate({heatmapdata()})
    results.data <- isolate({heatmapresultdata()})

    colmdata <- isolate({input$colmetadataheatmap})
    rowmdata <- isolate({input$rowmetadataheatmap})
    cluster_rows <- isolate({input$heatmap.clusterrows})
    show_column_names <- isolate({input$heatmap.show.col.names})
    show_row_names <- isolate({input$heatmap.show.row.names})
    cluster_columns <- isolate({input$heatmap.clustercol})
    sortCol <- isolate({input$heatmapSortCol})
    scale <- isolate({input$heatmapScale})
    title <-  isolate({input$heatmapMain})
    values.label <-  isolate({input$heatmapLabel})
    rownames.size <-  isolate({input$heatmapRownamesSize})

    if(isolate({input$heatmap.colorsCb})) {
        color.levels <- c(isolate({input$heatmapcolMin}),
                          isolate({input$heatmapcolMid}),
                          isolate({input$heatmapcolMax}))
    } else {
        color.levels <- NULL
    }
    if(isolate({input$heatmap.extremesCb})) {
        extrems <- c(isolate({input$heatmapExtremeMin}),
                     isolate({input$heatmapExtremeMid}),
                     isolate({input$heatmapExtremeMax}))
    } else {
        extrems <- NULL
    }
    if(isolate({input$heatmapTypeInputRb}) == "met") type <-  "methylation"
    if(isolate({input$heatmapTypeInputRb}) == "exp") type <-  "expression"

    if(nchar(sortCol) == 0 &  isolate({input$heatmap.sortCb})){
        createAlert(session, "heatmapmessage", "heatmapAlert", title = "Columns metadata", style =  "danger",
                    content = paste0("Please select the heatmapSortCol"),append = FALSE)
        return(NULL)
    }

    if(isolate({input$heatmapTypeInputRb})=="met"){
        # ---------------- probes selection
        if(isolate({input$heatmapProbesInputRb}) == "Status"){
            sig.probes <- ""
            if(isolate({input$heatmap.hypoprobesCb})) sig.probes <- c("Hypomethylated")
            if(isolate({input$heatmap.hyperprobesCb})) sig.probes <- c("Hypermethylated",sig.probes)
            sig.probes <- paste(sig.probes,"in",group2)
            # Get hypo methylated and hypermethylated probes
            group1.col <- gsub("[[:punct:]]| ", ".", group1)
            group2.col <- gsub("[[:punct:]]| ", ".", group2)
            idx <- paste("status",group1.col,group2.col, sep=".")
            probes <- apply(sapply(sig.probes, function(x) {grepl(x,results.data[,idx])}),1,any)

            if(length(probes) == 0){
                createAlert(session, "heatmapmessage", "heatmapAlert", title = "No significant probes", style =  "danger",
                            content = paste0("There are no significant probes"),append = FALSE)
                return(NULL)
            }
            data <- data[probes,]
            results.data <- results.data[probes,]
        } else if(isolate({input$heatmaprobes}) == "Selection"){
            data <- data[rownames(data) %in% isolate({input$heatmaprobes}) ]
        } else {
            sig.probes <- parse.textarea.input(isolate({input$heatmapProbesTextArea}))
            data <- data[rownames(data) %in% sig.probes,]
        }
    } else {
        if(isolate({input$heatmapGenesInputRb}) == "Status"){
            sig.genes <- ""
            if(isolate({input$heatmap.upGenesCb})) sig.genes <- c("Upregulated")
            if(isolate({input$heatmap.downGenewsCb})) sig.genes <- c("Downregulated",sig.genes)
            sig.genes <- paste(sig.genes,"in",group2)
            # Get hypo methylated and hypermethylated probes
            genes <-  gsub("[[:punct:]]| ", ".",results.data[,"status"]) %in%  gsub("[[:punct:]]| ", ".",sig.genes)
            results.data <- results.data[genes,]
            data <- data[rownames(data) %in% results.data$mRNA |rownames(data) %in% results.data$Gene_symbol ,]
        } else if(isolate({input$heatmapGenesInputRb}) == "Selection"){
            data <- data[rownames(data) %in% isolate({input$heatmagenes}) ]
        } else {
            sig.genes <- parse.textarea.input(isolate({input$heatmapGenesTextArea}))
            data <- data[rownames(data) %in% sig.genes]
        }

    }
    # ---------------- col.metadata
    if(!("barcode" %in% colnames(colData(data)))){
        createAlert(session, "heatmapmessage", "heatmapAlert", title = "Columns metadata", style =  "danger",
                    content = paste0("Sorry, but I need a barcode column to map the Summarized Experiment object",append = FALSE))
        return(NULL)
    }

    col.metadata <- NULL
    if(!is.null(colmdata)) {
        if(length(colmdata) > 0) col.metadata <- subset(colData(data), select=c("barcode",colmdata))
    }

    # ---------------- row.metadata
    row.metadata <- NULL
    if(!is.null(rowmdata)) {
        # We will consider the values in the resutls object, as the main object may not have the results
        if(length(rowmdata) > 0) row.metadata <- subset(results.data, select=c(rowmdata))
    }
    if(isolate({input$heatmaplog2_plus_one})) {
        assay(data) <- log2(assay(data) + 1) # if there were 0 values the scale was giving bugs
    }

    withProgress(message = 'Creating plot',
                 detail = 'This may take a while...', value = 0, {
                     if(!isolate({input$heatmap.sortCb})) {
                         p <-  TCGAvisualize_Heatmap(data=assay(data),
                                                     col.metadata=col.metadata,
                                                     row.metadata=row.metadata,
                                                     title = title,
                                                     color.levels = color.levels,
                                                     rownames.size = rownames.size,
                                                     values.label = values.label,
                                                     filename = NULL,
                                                     extrems = extrems,
                                                     cluster_rows = cluster_rows,
                                                     show_column_names = show_column_names,
                                                     cluster_columns = cluster_columns,
                                                     show_row_names = show_row_names,
                                                     type = type,
                                                     scale = scale)
                     } else {
                         p <-  TCGAvisualize_Heatmap(data=assay(data),
                                                     col.metadata=col.metadata,
                                                     row.metadata=row.metadata,
                                                     title = title,
                                                     extrems = extrems,
                                                     color.levels = color.levels,
                                                     rownames.size = rownames.size,
                                                     values.label = values.label,
                                                     filename = NULL,
                                                     cluster_rows = cluster_rows,
                                                     show_column_names = show_column_names,
                                                     cluster_columns = cluster_columns,
                                                     show_row_names = show_row_names,
                                                     sortCol = sortCol,
                                                     type = type,
                                                     scale = scale)
                     }
                     incProgress(1/2)
                     p
                 })
})
observeEvent(input$heatmapPlotBt , {
    output$heatmap.plotting <- renderPlot({
        closeAlert(session,"heatmapAlert")
        # get information from file
        ComplexHeatmap::draw(heatmap.plot())
    })})

observeEvent(input$heatmapPlotBt , {
    updateCollapse(session, "collapseHeatmap", open = "Heatmap")
    output$heatmapPlot <- renderUI({
        plotOutput("heatmap.plotting", width = paste0(isolate({input$heatmapwidth}), "%"), height = isolate({input$heatmapheight}))
    })
})


output$saveheatmappicture <- downloadHandler(
    filename = function(){input$heatmapPlot.filename},
    content = function(file) {
        if(tools::file_ext(file) == "png") {
            grDevices::png(file, width = 10, height = 10,
                           res = 300, units = "in")
        } else if(tools::file_ext(file) == "pdf") {
            grDevices::pdf(file, width = 10, height = 10)
        } else if(tools::file_ext(file) == "svg") {
            grDevices::svg(file, width = 10, height = 10)
        }
        ComplexHeatmap::draw(heatmap.plot())
        dev.off()
    })

