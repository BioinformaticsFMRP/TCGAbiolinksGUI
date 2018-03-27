##----------------------------------------------------------------------
#                             Volcano plot
##----------------------------------------------------------------------

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Handling visibility
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
shinyjs::hide("volcanoNamesFill")

observeEvent(input$volcanoNames, {
    if(input$volcanoNames){
        shinyjs::show("volcanoNamesFill")
    } else {
        shinyjs::hide("volcanoNamesFill")
    }
})
observeEvent(input$volcanoInputRb, {
    if(input$volcanoInputRb == "met"){
        shinyjs::show("colHypomethylated")
        shinyjs::show("colHypermethylated")
        shinyjs::hide("colUpregulated")
        shinyjs::hide("colDownregulated")
        shinyjs::show("volcanoxcutMet")
        shinyjs::hide("volcanoxcutExp")
    } else  if(input$volcanoInputRb == "exp"){
        shinyjs::hide("colHypomethylated")
        shinyjs::hide("colHypermethylated")
        shinyjs::show("colUpregulated")
        shinyjs::show("colDownregulated")
        shinyjs::show("volcanoxcutExp")
        shinyjs::hide("volcanoxcutMet")
    }
})
# If the user wants to hilght as probe or gene
observe({
    if(!is.null(input$volcanoHighlight)){
        updateCheckboxInput(session, "volcanoNames",  value = TRUE)
        updateCheckboxInput(session, "volcanoNamesFill",  value = TRUE)
        updateSelectizeInput(session, 'volcanoShowHighlitgh', selected = "highlighted")
    } else {
        updateSelectizeInput(session, 'volcanoShowHighlitgh', selected = "significant")
    }
})
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# File selection
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observe({
    shinyFileChoose(input, 'volcanofile', roots=get.volumes(input$workingDir),
                    session=session, restrictions=system.file(package='base'),
                    filetypes=c('excel', 'csv'))
})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# DATA INPUT
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
volcanodata <-  reactive({
    inFile <- input$volcanofile
    if (is.null(inFile)) return(NULL)
    file  <- as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), inFile)$datapath)
    # verify if the file is a csv
    ext <- tools::file_ext(file)
    if(ext != "csv"){
        createAlert(session, "dmrmessage", "dmrAlert", title = "Data input error", style = "danger",
                    content = paste0("Sorry, but I'm expecting a csv file, but I got a: ",
                                     ext), append = FALSE)
        return(NULL)
    }

    withProgress(message = 'Loading data',
                 detail = 'This may take a while...', value = 0, {
                     df <- as.data.frame(read_csv(file)); df$X1 <- NULL
                     incProgress(1, detail = "Completed")
                 })
    return(df)
})


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# UPDATING FIELDS AFTER DATA INPUT
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observeEvent(input$volcanofile, {
    file  <- basename(as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), input$volcanofile)$datapath))
    if(length(file) > 0){
        file <- unlist(str_split(file,"_"))
        group1 <- file[4]
        group2 <- file[5]
        pcut <- file[7]
        meancut <- gsub(".csv","",file[9])
        updateNumericInput(session, "volcanoxcutMet", value = meancut)
        updateNumericInput(session, "volcanoxcutExp", value = meancut)
        updateNumericInput(session, "volcanoycut", value = pcut)
    }
})
# automatically change the type based in the input
observe({
    if(!is.null(input$volcanofile)){
        file  <- basename(as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), input$volcanofile)$datapath))
        selected <- "met"
        if(grepl("DEA",file))  selected <- "exp"
        updateRadioButtons(session, "volcanoInputRb", selected = selected)
    }
})

# Update choises to select to highlight
observe({
    data <- volcanodata()
    if(!is.null(data)) {
        file  <- basename(as.character(parseFilePaths(get.volumes(isolate({input$workingDir})),
                                                      input$volcanofile)$datapath))
        if(grepl("DEA",file)){
            if("Gene_symbol" %in% colnames(data)){
                updateSelectizeInput(session, 'volcanoHighlight',
                                     choices = as.character(na.omit(unique(data$Gene_symbol))), server = TRUE)
            }
        } else {
            updateSelectizeInput(session, 'volcanoHighlight',
                                 choices = as.character(na.omit(unique(data$probeID))), server = TRUE)
        }
    }
})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# STATUS BOXES
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observeEvent(input$volcanoPlotBt , {
    output$volcanoBoxUp <- renderValueBox({
        ret <- isolate({volcano.values()})
        if(is.null(ret)) {
            value <- 0
        } else {
            value <- ret$up
        }
        valueBox(
            value = value,
            subtitle = ret$label[2],
            icon = icon("arrow-up"),
            color = "red"
        )
    })
})
observeEvent(input$volcanoPlotBt , {
    output$volcanoBoxInsig <- renderValueBox({
        ret <- isolate({volcano.values()})
        if(is.null(ret)) {
            value <- 0
        } else {
            value <- ret$insig
        }
        valueBox(
            value = value,
            ret$label[1],
            icon = icon("minus"),
            color = "black"
        )
    })
})
observeEvent(input$volcanoPlotBt , {
    output$volcanoBoxDown <- renderValueBox({
        ret <- isolate({volcano.values()})
        if(is.null(ret)) {
            value <- 0
        } else {
            value <- ret$down
        }
        valueBox(
            value = value,
            ret$label[3],
            icon = icon("arrow-down"),
            color = "olive"
        )
    })
})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# PLOT
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
volcano.values <- reactive({
    if(input$volcanoPlotBt){
        closeAlert(session, "volcanoAlert")

        # read csv file with results
        data <- volcanodata()
        if(is.null(data)) return(NULL)
        names.fill <- isolate({input$volcanoNamesFill})
        if(isolate({input$volcanoInputRb})=="met") {
            x.cut <- isolate({as.numeric(input$volcanoxcutMet)})
        } else {
            x.cut <- isolate({as.numeric(input$volcanoxcutExp)})
        }
        y.cut <- isolate({as.numeric(input$volcanoycut)})

        # Set parameters based in the filename
        # patterns are
        # DEA_result_groupCol_group1_group2_pcut_0.05_logFC.cut_0.csv
        # DMR_results_groupCol_group1_group2_pcut_0.05_meancut_0.3.csv
        file  <- basename(as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), input$volcanofile)$datapath))
        file <- unlist(str_split(file,"_"))
        groupCol <- file[3]
        group1 <- file[4]
        group2 <- file[5]
        names <- NULL

        # methylation pipeline
        if(isolate({input$volcanoInputRb})=="met"){
            group1.col <- gsub("[[:punct:]]| ", ".", group1)
            group2.col <- gsub("[[:punct:]]| ", ".", group2)
            diffcol <- paste("diffmean", group1.col, group2.col,sep = ".")
            pcol <- paste("p.value.adj", group1.col, group2.col,sep = ".")
            if(!(pcol %in% colnames(data) & diffcol %in% colnames(data) )) {
                createAlert(session, "volcanomessage", "volcanoAlert", title = "Error", style =  "success",
                            content = "We couldn't find the right columns in the data", append = FALSE)
            }
            if(isolate({input$volcanoNames})) names <- data$probeID
            label <- c("Not Significant",
                       "Hypermethylated",
                       "Hypomethylated")
            label[2:3] <-  paste(label[2:3], "in", group2)

            # Update data into a file

            statuscol <- paste("status",group1.col,group2.col,sep = ".")
            statuscol2 <- paste("status",group2.col,group1.col,sep = ".")
            data[,statuscol] <-  "Not Significant"
            data[,statuscol2] <-  "Not Significant"

            # get significant data
            sig <-  data[,pcol] < y.cut
            sig[is.na(sig)] <- FALSE
            # hypermethylated samples compared to old state
            hyper <- data[,diffcol]  > x.cut
            hyper[is.na(hyper)] <- FALSE

            # hypomethylated samples compared to old state
            hypo <- data[,diffcol] < (-x.cut)
            hypo[is.na(hypo)] <- FALSE

            if (any(hyper & sig)) data[hyper & sig,statuscol] <- paste("Hypermethylated","in", group2)
            if (any(hyper & sig)) data[hyper & sig,statuscol2] <- paste("Hypomethylated","in", group1)
            if (any(hypo & sig)) data[hypo & sig,statuscol] <- paste("Hypomethylated","in", group2)
            if (any(hypo & sig)) data[hypo & sig,statuscol2] <- paste("Hypermethylated","in", group1)
            insig.count <- nrow(data) - table(sig)["TRUE"]
            up.count <- table(hyper & sig)["TRUE"]
            down.count <- table(hypo & sig)["TRUE"]
            rownames(data) <- data$probeID
            if(isolate({input$volcanoSave})){
                getPath <- parseDirPath(get.volumes(isolate({input$workingDir})), input$workingDir)
                if (length(getPath) == 0) getPath <- paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI")
                csv <- paste0(paste("DMR_results",
                                    gsub("_",".",groupCol),
                                    group1.col,
                                    group2.col,
                                    "pcut",y.cut,
                                    "meancut",x.cut,
                                    sep = "_"),
                              ".csv")
                write_csv(data,path =  file.path(getPath, csv))
                createAlert(session, "volcanomessage", "volcanoAlert", title = "File created", style =  "success",
                            content = paste0(file.path(getPath, csv)), append = FALSE)

            }
            withProgress(message = 'Creating plot',
                         detail = 'This may take a while...', value = 0, {
                             p <-  TCGAVisualize_volcano(x = data[,diffcol],
                                                         y = data[,pcol],
                                                         ylab =   expression(paste(-Log[10],
                                                                                   " (FDR corrected -P values)")),
                                                         xlab =  expression(paste(
                                                             "DNA Methylation difference (",beta,"-values)")
                                                         ),
                                                         color = c(isolate({input$colinsignificant}),
                                                                   isolate({input$colHypermethylated}),
                                                                   isolate({input$colHypomethylated})),
                                                         title =  paste("Volcano plot", "(", group2, "vs", group1,")"),
                                                         legend=  "Legend",
                                                         label = label,
                                                         names = names,
                                                         names.fill = names.fill,
                                                         x.cut = x.cut,
                                                         y.cut = y.cut,
                                                         show.names = isolate({input$volcanoShowHighlitgh}),
                                                         highlight=isolate({input$volcanoHighlight}),
                                                         highlight.color = isolate({input$volcanoColHighlight}),
                                                         filename = NULL)
                         })
        } else {

            label <- c("Not Significant",
                       "Upregulated",
                       "Downregulated")
            label[2:3] <-  paste(label[2:3], "in", group2)
            if(isolate({input$volcanoNames})) names <- as.character(data$Gene_symbol)
            data$status <- "Insignificant"
            data[data$logFC >= x.cut & data$FDR <= y.cut,"status"] <- paste0("Upregulated in ", group2)
            data[data$logFC <= -x.cut & data$FDR <= y.cut,"status"] <- paste0("Downregulated in ", group2)

            up.count <- table(data$logFC >= x.cut & data$FDR <= y.cut)["TRUE"]
            if(is.na(up.count)) up.count <- 0
            down.count <- table(data$logFC <= -x.cut & data$FDR <= y.cut)["TRUE"]
            if(is.na(down.count)) down.count <- 0
            insig.count <-  nrow(data) -  down.count - up.count

            # Update data into a file
            if(isolate({input$volcanoSave})){
                getPath <- parseDirPath(get.volumes(isolate({input$workingDir})), input$workingDir)
                if (length(getPath) == 0) getPath <- paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI")

                out.filename <- paste0(paste("DEA_results",
                                             gsub("_",".",groupCol),
                                             gsub("_",".",group1),
                                             gsub("_",".",group2),
                                             "pcut", y.cut,
                                             "logFC.cut",x.cut,
                                             sep="_"),
                                       ".csv")
                out.filename <- file.path(getPath,out.filename)
                write_csv(data, path = out.filename)
                createAlert(session, "volcanomessage", "volcanoAlert", title = "File created", style =  "success",
                            content =  paste0(out.filename), append = FALSE)
            }

            withProgress(message = 'Creating plot',
                         detail = 'This may take a while...', value = 0, {
                             p <- TCGAVisualize_volcano(x = data$logFC,
                                                        y = data$FDR,
                                                        ylab =   expression(paste(-Log[10],
                                                                                  " (FDR corrected -P values)")),
                                                        xlab = " Gene expression fold change (Log2)",
                                                        color = c(isolate({input$colinsignificant}),
                                                                  isolate({input$colUpregulated}),
                                                                  isolate({input$colDownregulated})),
                                                        title =  paste("Volcano plot", "(", group2, "vs", group1,")"),
                                                        legend=  "Legend",
                                                        label = label,
                                                        names = names,
                                                        x.cut = x.cut,
                                                        y.cut = y.cut,
                                                        show.names = isolate({input$volcanoShowHighlitgh}),
                                                        highlight=isolate({input$volcanoHighlight}),
                                                        highlight.color = isolate({input$volcanoColHighlight}),
                                                        filename = NULL)
                         })
        }
    }
    ret <- list(plot = p, up = up.count, down = down.count, insig =insig.count, label = label)
})
observeEvent(input$volcanoPlotBt , {
    output$volcano.plot <- renderPlot({
        ret <- isolate({volcano.values()})
        if(is.null(ret)) return(NULL)
        ret$plot
    })
})
observeEvent(input$volcanoPlotBt , {
    updateCollapse(session, "collapseVolcano", open = "Volcano plot")
    output$volcanoPlot <- renderUI({
        plotOutput("volcano.plot", width = paste0(isolate({input$volcanowidth}), "%"), height = isolate({input$volcanoheight}))
    })})


output$savevolvanopicture <- downloadHandler(
    filename = function(){input$volcanoPlot.filename},
    content = function(file) {
        if(tools::file_ext(input$volcanoPlot.filename) == "png") {
            device <- function(..., width, height) {
                grDevices::png(..., width = 10, height = 10,
                               res = 300, units = "in")
            }
        } else if(tools::file_ext(input$volcanoPlot.filename) == "pdf") {
            device <- function(..., width, height) {
                grDevices::pdf(..., width = 10, height = 10)
            }
        } else if(tools::file_ext(input$volcanoPlot.filename) == "svg") {
            device <- function(..., width, height) {
                grDevices::svg(..., width = 10, height = 10)
            }
        }

        ggsave(file, plot = volcano.values()$plot, device = device)
    })
