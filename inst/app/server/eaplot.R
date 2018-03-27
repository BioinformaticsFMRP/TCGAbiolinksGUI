#--------------------------------------------------------------------
#                                    EA
#--------------------------------------------------------------------
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Handling visibility
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observeEvent(input$tcgaEaInputRb, {
    if(input$tcgaEaInputRb == "text") {
        shinyjs::show("eaGenesTextArea")
        shinyjs::hide("eagenes")
        shinyjs::hide("eaGenesFiles")
    } else if(input$tcgaEaInputRb == "Selection") {
        shinyjs::hide("eaGenesTextArea")
        shinyjs::show("eagenes")
        shinyjs::hide("eaGenesFiles")
    } else {
        shinyjs::hide("eaGenesTextArea")
        shinyjs::hide("eagenes")
        shinyjs::show("eaGenesFiles")
    }
})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# File selection
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observe({
    shinyFileChoose(input,
                    'eaGenesFiles',
                    roots=get.volumes(input$workingDir),
                    session=session,
                    restrictions=system.file(package='base'))
})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# DATA input
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
eaGenesByFile <- function(){
    inFile <- input$eaGenesFiles
    if (is.null(inFile)) return(NULL)
    file  <- as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), inFile)$datapath)
    if(tools::file_ext(file)=="csv"){
        df <- as.data.frame(read_csv(file));
        rownames(df) <- df[,1]
        df[,1] <- NULL
    } else if(tools::file_ext(file)=="rda"){
        df <- get(load(file))
    } else if(tools::file_ext(file)=="txt"){
        df <- read.table(file,header = T)
    } else {
        createAlert(session, "eamessage", "eaAlert", title = "Data input error", style =  "danger",
                    content = paste0("Sorry, but I'm expecting a csv, rda or txt file, but I got a: ",
                                     tools::file_ext(file)), append = FALSE)
        return(NULL)
    }
    genes <- NULL
    # if a data frame return the column with gene symbols
    if(class(df)==class(data.frame())){
        if("mRNA" %in% colnames(df)){
            df <- subset(df,df$status != "Insignificant")
            aux <- strsplit(df$mRNA,"\\|")
            genes <- unlist(lapply(aux,function(x) x[1]))
        } else if("Gene_symbol" %in% colnames(df)){
            genes <- df$Gene_symbol
        } else {
            createAlert(session, "eamessage", "eaAlert", title = "Data input error", style =  "danger",
                        content = paste0("Sorry, but I'm expecting a column called Gene_symbol "), append = FALSE)
            return(NULL)
        }
    }
    return(genes)
}
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# UPDATING FIELDS AFTER DATA INPUT
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observe({
    updateSelectizeInput(session,
                         'eagenes',
                         choices = unique(rownames(TCGAbiolinks:::EAGenes)),
                         server = TRUE)
})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# PLOT
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=

observeEvent(input$eaplot , {
    updateCollapse(session, "collapseEA", open = "EA plots")
    output$eaPlot <- renderUI({
        plotOutput("ea.plotting", width = paste0(isolate({input$eawidth}), "%"), height = isolate({input$eaheight}))
    })})

ea.plot  <-  reactive({
    input$eaplot
    textarea <- isolate({input$eaGenesTextArea})
    if(isolate({input$tcgaEaInputRb}) == "text" & !is.null(textarea)){
        genes <- toupper(parse.textarea.input(textarea))
        not.found <- genes[!(genes %in% TCGAbiolinks:::EAGenes$Gene)]
        if(length(not.found) > 0){
            closeAlert("eaAlert")
            createAlert(session, "eamessage", "eaAlert", title = "Data input error", style =  "danger",
                        content = paste0("Sorry, I cant't find these genes: ", not.found), append = FALSE)
            genes <-  genes[genes %in% TCGAbiolinks:::EAGenes$Gene]
        }
    } else if(isolate({input$tcgaEaInputRb}) == "Selection"){
        genes <- isolate({input$eagenes})
    } else{
        genes <- eaGenesByFile()
        not.found <- genes[!(genes %in% TCGAbiolinks:::EAGenes$Gene)]
        if(length(not.found) > 0){
            createAlert(session, "eamessage", "eaAlert", title = "Data input error", style =  "danger",
                        content = paste0("Sorry, I cant't find these genes: ", not.found), append = FALSE)
            genes <-  genes[genes %in% TCGAbiolinks:::EAGenes$Gene]
        }
    }

    xlim <- NULL
    if(isolate({input$eaxlim}) > 0) xlim <- c(0,isolate({input$eaxlim}))
    withProgress(message = 'Creating plot',
                 detail = 'This may take a while...', value = 0, {
                     ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes", genes)

                     ResMF <- NULL
                     ResBP <- NULL
                     ResCC <- NULL
                     ResPat <- NULL
                     if(length(grep("NA",ansEA$ResBP)) != ncol(ansEA$ResBP) &  isolate({input$eaPlotBPCB})) ResBP <- ansEA$ResBP
                     if(length(grep("NA",ansEA$ResCC)) != ncol(ansEA$ResCC) &  isolate({input$eaPlotCCCB})) ResCC <- ansEA$ResCC
                     if(length(grep("NA",ansEA$ResMF)) != ncol(ansEA$ResMF) &  isolate({input$eaPlotMFCB})) ResMF <- ansEA$ResMF
                     if(length(grep("NA",ansEA$ResPat)) != ncol(ansEA$ResPat) &  isolate({input$eaPlotPatCB})) ResPat <- ansEA$ResPat

                     plots <- table(c(is.null(ResMF),is.null(ResBP),is.null(ResCC),is.null(ResPat)))
                     if(!"FALSE" %in% names(plots) > 0){
                         createAlert(session, "eamessage", "eaAlert", title = "Data input error", style =  "danger",
                                     content = paste0("Sorry, no relevant results were found"), append = FALSE)
                     }
                     nbPlot <- plots["FALSE"]
                     if(nbPlot == 1)  mfrow = c(1,1)
                     if(nbPlot == 2)  mfrow = c(2,1)
                     if(nbPlot > 2)  mfrow = c(2,2)
                     if(is.na(nbPlot)) return (NULL)
                     # Enrichment Analysis EA (TCGAVisualize)
                     # Gene Ontology (GO) and Pathway enrichment barPlot
                     p <- tryCatch({
                         TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),
                                                 GOBPTab = ResBP,
                                                 GOCCTab = ResCC,
                                                 GOMFTab = ResMF,
                                                 PathTab = ResPat,
                                                 color = c(isolate({input$colBP}),
                                                           isolate({input$colCC}),
                                                           isolate({input$colMF}),
                                                           isolate({input$colPat})),
                                                 nRGTab = genes,
                                                 text.size = isolate({input$eaSizeText}),
                                                 xlim = xlim,
                                                 mfrow = mfrow,
                                                 nBar = isolate({input$nBar}),
                                                 filename = NULL)
                         p <- recordPlot(0)
                         return(p)
                     },  error = function(e) {
                         createAlert(session, "eamessage", "eaAlert", title = "No significant results", style =  "danger",
                                     content = paste0("The list of genes did not produce significant results"), append = FALSE)
                         return(NULL)
                     })
                 })
    return(p)
})

observeEvent(input$eaplot , {
    output$ea.plotting <- renderPlot({
        ea.plot()
    })
})


output$saveeapicture <- downloadHandler(
    filename = function(){input$eaPlot.filename},
    content = function(file) {
        if(tools::file_ext(file) == "png") {
            grDevices::png(file, width = 10, height = 10,
                           res = 300, units = "in")
        } else if(tools::file_ext(file) == "pdf") {
            grDevices::pdf(file, width = 10, height = 10)
        } else if(tools::file_ext(file) == "svg") {
            grDevices::svg(file, width = 10, height = 10)
        }
        print(ea.plot())
        dev.off()
    })


