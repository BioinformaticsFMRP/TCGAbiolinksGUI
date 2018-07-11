#----------------------------------------------------------------------
#   Oncoprint plot
#----------------------------------------------------------------------
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Handling visibility
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
shinyjs::hide("mafAnnotationcols")
shinyjs::hide("mafAnnotationpos")
shinyjs::hide("oncoGenes")
#shinyjs::hide("oncoInputRb")
shinyjs::hide("oncoGenesTextArea")

observeEvent(input$mafAnnotation, {
    if(!is.null(annotation.maf())){
        shinyjs::show("mafAnnotationcols")
        shinyjs::show("mafAnnotationpos")
    }
})

# The used added a file, now we can show the genes
observeEvent(input$maffile, {
    if(!is.null(mut.oncoprint())){
        shinyjs::show("oncoInputRb")
        shinyjs::show("oncoGenes")
    }
})
observeEvent(input$oncoInputRb, {
    if(input$oncoInputRb == "Selection"){
        shinyjs::hide("oncoGenesTextArea")
        shinyjs::hide("oncoGenesFiles")
        shinyjs::show("oncoGenes")
    } else  if(input$oncoInputRb == "text"){
        shinyjs::show("oncoGenesTextArea")
        shinyjs::hide("oncoGenes")
        shinyjs::hide("oncoGenesFiles")
    } else {
        shinyjs::hide("oncoGenesTextArea")
        shinyjs::hide("oncoGenes")
        shinyjs::show("oncoGenesFiles")
    }
})

observeEvent(input$oncoInformation, {
    if(input$oncoInformation == "Variant_Classification"){
        shinyjs::hide("colDEL")
        shinyjs::hide("colINS")
        shinyjs::hide("colSNP")
        shinyjs::hide("colDNP")
        shinyjs::show("colFSD")
        shinyjs::show("colFSI")
        shinyjs::show("colIFD")
        shinyjs::show("colIFI")
        shinyjs::show("colMM")
        shinyjs::show("colNM")
        shinyjs::show("colNSM")
        shinyjs::show("colRNA")
        shinyjs::show("colSIL")
        shinyjs::show("colSS")
        shinyjs::show("colTSS")
        shinyjs::show("colTR")
    } else {
        shinyjs::show("colDEL")
        shinyjs::show("colINS")
        shinyjs::show("colSNP")
        shinyjs::show("colDNP")
        shinyjs::hide("colFSD")
        shinyjs::hide("colFSI")
        shinyjs::hide("colIFD")
        shinyjs::hide("colIFI")
        shinyjs::hide("colMM")
        shinyjs::hide("colNM")
        shinyjs::hide("colNSM")
        shinyjs::hide("colRNA")
        shinyjs::hide("colSIL")
        shinyjs::hide("colSS")
        shinyjs::hide("colTSS")
        shinyjs::hide("colTR")
    }
})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# UPDATING FIELDS AFTER DATA INPUT
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# After user uploaded he MAF file we will add the choices availble
observe({
    updateSelectizeInput(session, 'oncoGenes', choices = unique(as.character(mut.oncoprint()$Hugo_Symbol)), server = TRUE)
})

observe({
    updateSelectizeInput(session, 'mafAnnotationcols',
                         choices = as.character(colnames(annotation.maf())), server = TRUE)
})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# File selection
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observe({
    shinyFileChoose(input, 'maffile', roots=get.volumes(input$workingDir), session=session,
                    restrictions=system.file(package='base'), filetypes=c('', "rda",'maf',"csv","maf.gz"))
    shinyFileChoose(input, 'mafAnnotation', roots=get.volumes(input$workingDir), session=session,
                    restrictions=system.file(package='base'), filetypes=c('', 'csv','rda'))
    shinyFileChoose(input, 'oncoGenesFiles', roots=get.volumes(input$workingDir), session=session, restrictions=system.file(package='base'))
})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# DATA INPUT
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
annotation.maf <- reactive({
    inFile <- input$mafAnnotation
    if(class(inFile) != "list") return(NULL)
    file  <- as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), input$mafAnnotation)$datapath)
    if(tools::file_ext(file)=="csv"){
        se <- as.data.frame(read_csv(file)); se$X1 <- NULL
    } else if(tools::file_ext(file)=="rda"){
        se <- get(load(file))
    }

    if(class(se)!= class(data.frame())){
        createAlert(session, "oncomessage", "oncoAlert", title = "Data input error", style =  "danger",
                    content = paste0("Sorry, but I'm expecting a Data frame object, but I got a: ",
                                     class(se)), append = FALSE)
        return(NULL)
    }
    return(se)
})
mut.oncoprint <-  reactive({
    inFile <- input$maffile
    if(class(inFile) != "list") return(NULL)
    inFile <- parseFilePaths(get.volumes(isolate({input$workingDir})), inFile)
    if (nrow(inFile) == 0) return(NULL)
    file <- as.character(inFile$datapath)
    withProgress(message = 'Reading MAF file',
                 detail = 'This may take a while...', value = 0, {
                     if(grepl("\\.maf\\.gz", file)){
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

genesByFile <- function(){
    inFile <- input$oncoGenesFiles
    if(class(inFile) != "list") return(NULL)
    file  <- as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), inFile)$datapath)
    if(tools::file_ext(file) == "csv"){
        se <- read_csv(file);
        rownames(df) <- df[,1]
        df[,1] <- NULL
    } else if(tools::file_ext(file) == "rda") {
        df <- get(load(file))
    } else if(tools::file_ext(file) == "txt") {
        df <- read.table(file,header = T)
    } else {
        createAlert(session, "oncomessage", "oncoAlert", title = "Data input error", style =  "danger",
                    content = paste0("Sorry, but I'm expecting a csv, rda or txt file, but I got a: ",
                                     tools::file_ext(file)), append = FALSE)
        return(NULL)
    }
    genes <- NULL
    # if a data frame return the column with gene symbols
    if(class(df) == class(data.frame())){
        if ("mRNA" %in% colnames(df)){
            if ("status" %in% colnames(df)) df <- subset(df,df$status != "Insignificant")
            aux <- strsplit(df$mRNA,"\\|")
            genes <- unlist(lapply(aux,function(x) x[1]))
        } else if ("Gene_symbol" %in% colnames(df)){
            genes <- df$Gene_symbol
        } else {
            createAlert(session, "oncomessage", "oncoAlert", title = "Data input error", style =  "danger",
                        content = paste0("Sorry, but I'm expecting a column called Gene_symbol "), append = FALSE)
            return(NULL)
        }
    }
    return(genes)
}


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Plot
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observeEvent(input$oncoprintPlot , {
    closeAlert(session,"oncoAlert")
    output$oncoploting <- renderPlot({
        mut <- isolate({mut.oncoprint()})
        annotation <- isolate({annotation.maf()})
        cols <- isolate({input$mafAnnotationcols})
        rm.empty.cols <- isolate({input$oncoRmCols})
        show.col.names <- isolate({input$oncoShowColsNames})
        textarea <- isolate({input$oncoGenesTextArea})
        show.row.barplot <- isolate({input$oncoShowRowBarplot})

        msg <- ""
        not.found <- c()

        if(isolate({input$oncoInputRb}) == "text" & !is.null(textarea)){
            genes <- toupper(parse.textarea.input(textarea))
            not.found <- genes[!(genes %in% mut$Hugo_Symbol)]
            if(length(not.found) > 0){
                msg <- paste0("Sorry, I cant't find these genes: ", paste(not.found,collapse = " "))
                genes <-  genes[genes %in% mut$Hugo_Symbol]
            }
        } else if(isolate({input$oncoInputRb}) == "Selection") {
            genes <- isolate({input$oncoGenes})
        } else if(isolate({input$oncoInputRb}) == "file") {
            genes <- genesByFile()
            not.found <- genes[!(genes %in% mut$Hugo_Symbol)]
            if(length(not.found) > 0){
                msg <- paste0("Sorry, I cant't find these genes: ", paste(not.found,collapse = " "))
                genes <-  genes[genes %in% mut$Hugo_Symbol]
            }
        }

        if(is.null(genes)){
            createAlert(session, "oncomessage", "oncoAlert", title = "Error", style =  "danger",
                        content = "Please select the genes (max 50)", append = TRUE)
            return(NULL)
        } else if( length(genes) > 200){
            createAlert(session, "oncomessage", "oncoAlert", title = "Errors", style =  "danger",
                        content = paste0("The limit of the genes is 200 \n You gave me: ",length(genes)), append = TRUE)
            return(NULL)
        } else if(all(genes == "")){
            createAlert(session, "oncomessage", "oncoAlert", title = "Errors", style =  "danger",
                        content = "From the list given none are mutated", append = TRUE)
            return(NULL)

        } else if(length(not.found) > 0){
            createAlert(session, "oncomessage", "oncoAlert", title = "Errors", style =  "danger",
                        content = msg, append = TRUE)
        }

        if(is.null(cols)) {
            annotation <- NULL
        } else if("bcr_patient_barcode" %in% colnames(annotation)) {
            annotation <- annotation[,c("bcr_patient_barcode",cols)]
        } else if("patient" %in% colnames(annotation)) {
            annotation <- annotation[,c("patient",cols)]
            colnames(annotation)[which(colnames(annotation) == "patient")] <- "bcr_patient_barcode"
        } else {
            createAlert(session, "oncomessage", "oncoAlert", title = "Error", style =  "danger",
                        content = "I couldn't find the bcr_patient_barcode or patient column in the annotation", append = TRUE)
            return(NULL)
        }

        if(is.null(mut)){
            createAlert(session, "oncomessage", "oncoAlert", title = "Error", style =  "danger",
                        content = "Please select a file", append = TRUE)
            return(NULL)
        }

        if(!is.null(annotation)){
            if(!"bcr_patient_barcode" %in% colnames(annotation)){
                createAlert(session, "oncomessage", "oncoAlert", title = "Error", style =  "danger",
                            content = "bcr_patient_barcode column should be in the annotation", append = TRUE)
                return(NULL)
            }
            idx <- match(substr(mut$Tumor_Sample_Barcode,1,12),annotation$bcr_patient_barcode)
            if(all(is.na(idx))){
                createAlert(session, "oncomessage", "oncoAlert", title = "Error", style =  "danger",
                            content = "Annotation samples does not match mutation", append = TRUE)
                return(NULL)
            }
        }

        withProgress(message = 'Creating plot',
                     detail = 'This may take a while...', value = 0, {

                         if(isolate({input$oncoInformation}) == "Variant_Type") {
                             color = c("background" = isolate(input$colBG),
                                       "SNP"=isolate(input$colSNP),"INS"=isolate(input$colINS),
                                       "DEL"=isolate(input$colDEL),"DNP"=isolate(input$colDNP))
                         } else {
                             color = c("background" = isolate(input$colBG),
                                       "Frame_Shift_Del"=isolate(input$colFSD),
                                       "Frame_Shift_Ins"=isolate(input$colFSI),
                                       "In_Frame_Del"=isolate(input$colIFD),
                                       "In_Frame_Ins"=isolate(input$colIFI),
                                       "Missense_Mutation"=isolate(input$colMM),
                                       "Nonsense_Mutation"=isolate(input$colNM),
                                       "Nonstop_Mutation"=isolate(input$colNSM),
                                       "RNA"=isolate(input$colRNA),
                                       "Silent"=isolate(input$colSIL),
                                       "In_Frame_Ins"=isolate(input$colIFI),
                                       "Splice_Site"=isolate(input$colSS),
                                       "Translation_Start_Site"=isolate(input$colTSS),
                                       "Targeted_Region"=isolate(input$colTR))
                         }
                         TCGAvisualize_oncoprint(mut=mut,genes=genes,annotation = annotation,
                                                 annotation.position=isolate(input$mafAnnotationpos),
                                                 rm.empty.columns = rm.empty.cols,
                                                 information = isolate({input$oncoInformation}),
                                                 show.column.names = show.col.names,
                                                 show.row.barplot = show.row.barplot,
                                                 label.font.size = isolate(input$oncoTextLabelSize),
                                                 rows.font.size = isolate(input$oncoTextRowSize),
                                                 dist.row =  isolate(input$oncoHSpace),
                                                 dist.col =  isolate(input$oncoWSpace),
                                                 row.order =  isolate(input$oncoRowSort),
                                                 annotation.legend.side = isolate(input$oncoAnnotationLegendSide),
                                                 heatmap.legend.side = isolate(input$oncoHeatmapLegendSide),
                                                 color = color)

                     })
    })
})

observeEvent(input$oncoprintPlot , {
    updateCollapse(session, "collapseOnco", open = "Oncoprint")
    output$oncoPlot <- renderUI({
        plotOutput("oncoploting", width = paste0(isolate({input$oncowidth}), "%"),
                   height = isolate({input$oncoheight}))
    })
})
