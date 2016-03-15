library(shiny)
library(shinyFiles)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(shinyBS)
options(shiny.maxRequestSize=300*1024^2)


#' @title  Server side
#' @description Server side - Download data from roadmap project
#' @param input - input signal
#' @param output - output signal
#' @importFrom downloader download
#' @keywords internal
biOMICsServer <- function(input, output, session) {
    setwd(Sys.getenv("HOME"))
    volumes <- c('Working directory'=getwd())
    shinyDirChoose(input, 'folder', roots=volumes, session=session, restrictions=system.file(package='base'))
    shinyDirChoose(input, 'reportfolder', roots=volumes, session=session, restrictions=system.file(package='base'))
    output$reportdirectorypath <- renderText({parseDirPath(volumes, input$reportfolder)})

    observeEvent(input$folder , {
        closeAlert(session, "ontdownloaddirAlert")
        createAlert(session, "ontdownloaddirmessage", "ontdownloaddirAlert", title = "Download folder", style =  "success",
                    content =   parseDirPath(volumes, input$folder), append = TRUE)
    })
    observeEvent(input$reportfolder , {
        closeAlert(session, "ontreportdirAlert")
        createAlert(session, "ontdownloaddirmessage", "ontreportdirAlert", title = "Report folder", style =  "success",
                    content =  parseDirPath(volumes, input$reportfolder), append = TRUE)
    })

    #----------------- Ontology
    observeEvent(input$ontSearchDownloadBt , {
        selected.rows <- isolate({input$allRows})
        if( length(selected.rows) == 0) {
            createAlert(session, "alert", "exampleAlert", title = "Invalid selection", style =  "danger",
                        content = "Please select the samples to download.", append = FALSE)
            return()
        } else {
            closeAlert(session, "exampleAlert")
        }

        getPath <- parseDirPath(volumes, isolate({input$folder}))
        getFtype <- input$ontftypeFilter

        link <- c()

        df <- data.frame(database =selected.rows[seq(1, length(selected.rows) , 5)],
                         ID =selected.rows[seq(2, length(selected.rows) , 5)],
                         Sample = selected.rows[seq(3, length(selected.rows) , 5)],
                         Experiment = selected.rows[seq(4, length(selected.rows) , 5)],
                         organism = selected.rows[seq(5, length(selected.rows) , 5)])
        withProgress(message = 'Downloading data',
                     detail = 'This may take a while...', value = 0, {
                         if(length(getPath) > 0) {
                             biOmicsDownload(df, path = getPath, enc.file.type = getFtype, rmap.file.type = getFtype)
                         } else {
                             biOmicsDownload(df, enc.file.type = getFtype, rmap.file.type = getFtype)
                         }

                     })
        createAlert(session, "alert", "exampleAlert", title = "Download completed", style =  "info",
                    content = "Your download has been completed.", append = FALSE)
    })

    dataInput <- reactive({
        withProgress(message = 'Searching in progress',
                     detail = 'This may take a while...', value = 0, {
                         indexes <- c()
                         query <- data.frame()
                         link <- c()
                         term <- isolate({input$ontSamplesFilter})
                         exp <-  isolate({input$ontExpFilter})
                         query <- biOmicsSearch(term, experiment = exp)
                         # improve using subset - subset(data,selection,projection)
                         biosample.encode <- get.obj("biosample.encode")
                         return(list(system = names(sort(table(biosample.encode[biosample.encode$biosample %in% query$Sample,]$system),decreasing = T)[1]),
                                     result = query))
                     })})

    observeEvent(input$ontReport, {
        result <- dataInput()$result

        print(result)
        if(is.null(result)) {
            createAlert(session, "alert", "exampleAlert", title = "No system", style =  "danger",
                        content = "No system was found. I can't create a report.", append = FALSE)
            return()

        }
        else if(is.null(input$reportfolder)) {
            createAlert(session, "alert", "exampleAlert", title = "No folder", style =  "danger",
                        content = "Please select a folder to create the report", append = FALSE)
            return()
        }  else {
            closeAlert(session, "exampleAlert")
        }

        withProgress(message = 'Creating report...',
                     detail = 'This may take a while...', value = 0, {
                         #create.report(dataInput()$result, system = dataInput()$system)
                         create.report(result, path = file.path(parseDirPath(volumes, input$reportfolder),"report"), system = dataInput()$system)
                     })
    })

    observeEvent(input$tcgaPrepareBt, {

        # read the data from the downloaded path
        # prepare it

        # Files types
        ftype <- NULL
        rnaseqFtype <- isolate({input$tcgaFrnaseqtypeFilter})
        rnaseqv2Ftype <- isolate({input$tcgaFrnaseqv2typeFilter})
        gwsFtype <- isolate({input$tcgaFgwstypeFilter})

        # Dir to saved the files
        getPath <- parseDirPath(volumes, input$tcgafolder)
        if (length(getPath) == 0) getPath <- "."
        samplesType <- input$tcgasamplestypeFilter

        save.dir <- parseDirPath(volumes, input$tcgapreparefolder)
        if(length(save.dir) == 0) {
            filename <- isolate({input$tcgafilename})
        } else {
            filename <- file.path(save.dir,isolate({input$tcgafilename}))
        }

        withProgress( message = 'Prepare in progress',
                      detail = 'This may take a while...', value = 0, {
                          df <- data.frame(name = input$allRows[seq(6, length(input$allRows), 7)])
                          x <- TCGAquery()
                          x <- x[x$name %in% df$name,]
                          if(length(unique(x$Platform)) > 1 & !all(grepl("humanmethylation",unique(x$Platform),ignore.case = TRUE))) {
                              print("We can't prepare these data together")
                              return(NULL)
                          }
                          ftype <- NULL
                          if ("IlluminaHiSeq_RNASeqV2" %in% unique(x$Platform)) ftype <- rnaseqv2Ftype
                          if ("IlluminaHiSeq_RNASeq"  %in% unique(x$Platform)) ftype <- rnaseqFtype
                          if ("Genome_Wide_SNP_6" %in% unique(x$Platform)) ftype <- gwsFtype
                          if (length(samplesType) == 0) {
                              samples <- NULL
                          } else {
                              samples <- unlist(lapply(samplesType,function(type){
                                  s <- unlist(str_split(x$barcode,","))
                                  s[grep(type,substr(s,14,15))]
                              }))
                          }

                          trash <- TCGAprepare(x,dir = getPath,
                                               summarizedExperiment = isolate({as.logical(input$prepareRb)}),
                                               save = TRUE,
                                               type = ftype,
                                               filename=filename,
                                               samples = samples,
                                               add.subtype = isolate({input$addSubTypeTCGA}))
                      })
        closeAlert(session, "tcgaAlert")
        createAlert(session, "tcgasearchmessage", "tcgaAlert", title = "Prepare completed", style =  "success",
                    content =  paste0("Saved in: ", getPath), append = FALSE)
    })

    observeEvent(input$ontSearchBt, {

        if(nchar(input$ontSamplesFilter) < 3) {
            closeAlert(session, "exampleAlert")
            createAlert(session, "alert", "exampleAlert", title = "Invalid term", style =  "danger",
                        content = "Input should be creater than 3 characters.", append = FALSE)
            return()
        } else {
            closeAlert(session, "exampleAlert")
        }
        output$ontSearchtbl <- renderDataTable({
            data <- isolate({dataInput()})
            createAlert(session, "alert", "exampleAlert", title = "Term mapped to the system:", style =  "info",
                        content = data$system , append = FALSE)
            result <- data$result
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
                           "sSwfPath" = paste0("//cdnjs.cloudflare.com/ajax/",
                                               "libs/datatables-tabletools/",
                                               "2.2.3/swf/copy_csv_xls.swf"),
                           "aButtons" = list(
                               list("sExtends" = "collection",
                                    "sButtonText" = "Save",
                                    "aButtons" = c("csv","xls")
                               )
                           )
                       )
        ), callback = "function(table) {table.on('click.dt', 'tr', function() {Shiny.onInputChange('allRows',table.rows('.selected').data().toArray());});}"
        )})

    #------------- TCGA ------------------

    # Table render
    observeEvent(input$tcgaSearchBt, {
        output$tcgaSearchtbl <- renderDataTable({
            tumor = input$tcgaTumorFilter
            platform = input$tcgaExpFilter
            level = input$tcgaLevelFilter

            tbl <- data.frame()
            if(length(level) > 0){
                for(i in level){
                    tbl <- rbind(tbl,
                                 TCGAquery(tumor = tumor,
                                           platform = platform,
                                           level = i))
                }
            } else {
                tbl <- TCGAquery(tumor = tumor,
                                 platform = platform,
                                 level = level)
            }
            if(is.null(tbl)){
                createAlert(session, "tcgasearchmessage", "tcgasearchAlert", title = "No results found", style =  "warning",
                            content = "Sorry there are no results for your query.", append = FALSE)
                return()
            } else if(nrow(tbl) ==0) {
                createAlert(session, "tcgasearchmessage", "tcgasearchAlert", title = "No results found", style =  "warning",
                            content = "Sorry there are no results for your query.", append = FALSE)
                return()
            } else {
                closeAlert(session, "tcgasearchAlert")
                tbl$Level <- stringr::str_extract(tbl$name,"Level_[1-3]|mage-tab|aux")
                tbl <- tbl[,c(1,10,11,12,13,7,8)]
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
        )})

    # Download
    shinyDirChoose(input, 'tcgafolder', roots=volumes, session=session,
                   restrictions=system.file(package='base'))
    shinyDirChoose(input, 'tcgapreparefolder', roots=volumes, session=session,
                   restrictions=system.file(package='base'))

    observeEvent(input$tcgafolder , {
        closeAlert(session, "tcgadownloaddirAlert")
        createAlert(session, "tcgaddirmessage", "tcgadownloaddirAlert", title = "Download folder", style =  "success",
                    content =  parseDirPath(volumes, input$tcgafolder), append = TRUE)
    })

    observeEvent(input$tcgapreparefolder , {
        closeAlert(session, "tcgapreparedirAlert")
        createAlert(session, "tcgaddirmessage", "tcgapreparedirAlert", title = "Prepare folder", style =  "success",
                    content =  parseDirPath(volumes, input$tcgapreparefolder), append = TRUE)
    })


    observeEvent(input$tcgaDownloadBt,{

        # Files types
        ftype <- NULL
        rnaseqFtype <- isolate({input$tcgaFrnaseqtypeFilter})
        rnaseqv2Ftype <- isolate({input$tcgaFrnaseqv2typeFilter})
        gwsFtype <- isolate({input$tcgaFgwstypeFilter})

        # Dir to save the files
        getPath <- parseDirPath(volumes, input$tcgafolder)
        if (length(getPath) == 0) getPath <- "."
        samplesType <- input$tcgasamplestypeFilter


        withProgress(message = 'Download in progress',
                     detail = 'This may take a while...', value = 0, {
                         df <- data.frame(name = input$allRows[seq(6, length(input$allRows), 7)])
                         x <- TCGAquery()
                         x <- x[x$name %in% df$name,]

                         for (i in 1:nrow(x)) {
                             incProgress(1/nrow(x))
                             aux <- x[i,]
                             if (aux$Platform == "IlluminaHiSeq_RNASeqV2") ftype <- rnaseqv2Ftype
                             if (aux$Platform == "IlluminaHiSeq_RNASeq") ftype <- rnaseqFtype
                             if (aux$Platform == "Genome_Wide_SNP_6") ftype <- gwsFtype
                             if (length(samplesType) == 0) {
                                 samples <- NULL
                             } else {
                                 samples <- unlist(lapply(samplesType,function(type){
                                     s <- unlist(str_split(aux$barcode,","))
                                     s[grep(type,substr(s,14,15))]
                                 }))
                             }
                             TCGAdownload(x, path = getPath,type = ftype,samples = samples)
                         }})
        createAlert(session, "tcgasearchmessage", "tcgaAlert", title = "Download completed", style =  "success",
                    content =  paste0("Saved in: ", getPath), append = FALSE)
    })

    # Subtype

    observeEvent(input$tcgaSubtypeBt, {
        output$tcgaSearchtbl <- renderDataTable({
            tumor <- isolate({input$tcgasubtypeFilter})
            tbl <- data.frame()


            result = tryCatch({
                tbl <- rbind(tbl, TCGAquery_subtype(tumor = tumor))


            }, error = function(e) {
                createAlert(session, "tcgasearchmessage", "tcgasearchAlert", title = "No results found", style =  "warning",
                            content = "Sorry there are subtypes for your query.", append = FALSE)
            })

            if(is.null(tbl)){
                createAlert(session, "tcgasearchmessage", "tcgasearchAlert", title = "No results found", style =  "warning",
                            content = "Sorry there are subtypes for your query.", append = FALSE)
                return()
            } else if(nrow(tbl) ==0) {
                createAlert(session, "tcgasearchmessage", "tcgasearchAlert", title = "No results found", style =  "warning",
                            content = "Sorry there are subtypes for your query.", append = FALSE)
                return()
            } else {
                closeAlert(session, "tcgasearchAlert")
                doi <- c("aml"="doi:10.1056/NEJMoa1301689",
                         "blca"="doi:10.1038/nature12965",
                         "brca"="doi:10.1038/nature11412",
                         "coad"="doi:10.1038/nature11252",
                         "gbm"="doi:10.1016/j.cell.2015.12.028",
                         "lgg"="doi:10.1016/j.cell.2015.12.028",
                         "hnsc"="doi:10.1038/nature14129",
                         "kich"="doi:10.1016/j.ccr.2014.07.014",
                         "kirc"="doi:10.1038/nature12222",
                         "kirp"="doi:10.1056/NEJMoa1505917",
                         "lihc"="",
                         "luad"="doi:10.1038/nature13385",
                         "lusc"="doi:10.1038/nature11404",
                         "ovca"= "doi:10.1038/nature10166",
                         "pancan"="doi:10.1016/j.cell.2014.06.049",
                         "prad"="doi:10.1016/j.cell.2015.10.025",
                         "skcm"="doi:10.1016/j.cell.2015.05.044",
                         "stad"="doi:10.1038/nature13480",
                         "thca"="doi:10.1016/j.cell.2014.09.050",
                         "ucec"="doi:10.1038/nature12113",
                         "ucs"="")
                createAlert(session, "tcgasearchmessage", "tcgasearchAlert", title = "Source of the data", style =  "warning",
                            content = paste0(doi[tumor]), append = FALSE)
                if (isolate({input$saveSubtype})) save(tbl, file = paste0(tumor,"_subtype.rda"))
                return(tbl)
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
        )})

    observeEvent(input$tcgaClinicalBt, {
        output$tcgaSearchtbl <- renderDataTable({
            tumor <- isolate({input$tcgatumorClinicalFilter})
            type <- isolate({input$tcgaClinicalFilter})

            tbl <- data.frame()
            withProgress( message = 'Prepare in progress',
                          detail = 'This may take a while...', value = 0, {
                              result = tryCatch({
                                  tbl <- rbind(tbl, TCGAquery_clinic(tumor = tumor,clinical_data_type = type))

                              }, error = function(e) {
                                  createAlert(session, "tcgasearchmessage", "tcgasearchAlert", title = "No results found", style =  "warning",
                                              content = "Sorry there are clinical files for your query.", append = FALSE)
                              })
                          })
            if(is.null(tbl)){
                createAlert(session, "tcgasearchmessage", "tcgasearchAlert", title = "No results found", style =  "warning",
                            content = "Sorry there are clinical files for your query.", append = FALSE)
                return()
            } else if(nrow(tbl) ==0) {
                createAlert(session, "tcgasearchmessage", "tcgasearchAlert", title = "No results found", style =  "warning",
                            content = "Sorry there are clinical files for your query.", append = FALSE)
                return()
            } else {
                closeAlert(session, "tcgasearchAlert")
                if (isolate({input$saveClinical})) save(tbl, file = paste0(tumor,"_clinic_",type,".rda"))
                return(tbl)
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
        )})

    #------------- MAF

    # Table render
    output$maftbl <- renderDataTable({
        maf.files <- get.obj("maf.files")
        maf.files[,c(10,1,2,4,5,7)]
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

    # Download
    shinyDirChoose(input, 'maffolder', roots=volumes, session=session,
                   restrictions=system.file(package='base'))

    observeEvent(input$maffolder , {
        closeAlert(session, "oncoddirAlert")
        createAlert(session, "oncoddirmessage", "oncoddirAlert", title = "Download folder", style = "success",
                    content =  parseDirPath(volumes, input$maffolder), append = TRUE)
    })

    observeEvent(input$mafDownloadBt , {
        maf.files <- get.obj("maf.files")
        # Dir to save the files
        getPath <- parseDirPath(volumes, input$maffolder)
        if (length(getPath) == 0) getPath <- "."

        withProgress(message = 'Download in progress',
                     detail = 'This may take a while...', value = 0, {
                         df <- data.frame(name = input$allRows[seq(5, length(input$allRows), 7)])
                         print(df)
                         df <-   maf.files[maf.files$Archive.Name %in% df$name,]

                         for (i in 1:nrow(df)) {
                             incProgress(1/nrow(df))
                             fout <- file.path(getPath,basename(df[1,]$Deploy.Location))
                             if (!file.exists(fout))  downloader::download(df[1,]$Deploy.Location,fout)
                         }})
        closeAlert(session, "oncoAlert")
        createAlert(session, "oncomessage", "oncoAlert", title = "Download completed", style = "success",
                    content =  paste0("Saved file: ",fout), append = TRUE)
    })
    shinyFileChoose(input, 'maffile', roots=volumes, session=session, restrictions=system.file(package='base'))

    mut <- function(){
        #inFile <- input$maffile
        #if (is.null(inFile)) return(NULL)
        inFile <- parseFilePaths(volumes, input$maffile)
        if (nrow(inFile) == 0) return(NULL)
        file <- as.character(inFile$datapath)
        ret <- read.table(file, fill = TRUE,
                          comment.char = "#", header = TRUE, sep = "\t", quote='')
        return(ret)

    }
    observeEvent(input$oncoprintPlot , {
        output$oncoploting <- renderPlot({
            mut <- isolate({mut()})
            if(is.null(mut)){
                createAlert(session, "oncomessage", "oncoAlert", title = "Error", style =  "danger",
                            content = "Please select a file", append = TRUE)
                return(NULL)
            } else{
                closeAlert(session, "oncoAlert")
            }
            withProgress(message = 'Creating plot',
                         detail = 'This may take a while...', value = 0, {

                             create.oncoprint(mut=mut,genes=isolate(input$oncoGenes),
                                              color = c("SNP"=isolate(input$colSNP),"INS"=isolate(input$colINS),
                                                        "DEL"=isolate(input$colDEL),"DNP"=isolate(input$colDNP)))

                         })
        })})
    observeEvent(input$oncoprintPlot , {
        updateCollapse(session, "collapseOnco", open = "Oncoprint")
        output$oncoPlot <- renderUI({
            plotOutput("oncoploting", width = paste0(isolate({input$oncowidth}), "%"), height = isolate({input$oncoheight}))
        })})

    observe({
        updateSelectizeInput(session, 'oncoGenes', choices = as.character(mut()$Hugo_Symbol), server = TRUE)
    })
    observe({
        updateSelectizeInput(session, 'ontExpFilter', choices = as.character(get.obj("platforms")$Standard), server = TRUE)
    })
    observe({
        updateSelectizeInput(session, 'ontSamplesFilter', choices = as.character(c("",
                                                                                   union(
                                                                                       union(
                                                                                           get.obj("roadmap.db")$Sample.Name,
                                                                                           get.obj("encode.db")$biosample
                                                                                       ),
                                                                                       TCGAbiolinks::TCGAquery()$Disease
                                                                                   ))), server = TRUE)
    })
    observe({
        updateSelectizeInput(session, 'ontftypeFilter', choices = as.character(unique(get.obj("encode.db.files")$file_format)), server = TRUE)
    })

    ## DMR analysis
    observeEvent(input$dmrAnalysis , {

        # read the data from the downloaded path
        # prepare it
        se <- isolate({dmrdata()})
        se <- subset(se,subset = (rowSums(is.na(assay(se))) == 0))
        withProgress(message = 'DMR analysis in progress',
                     detail = 'This may take a while...', value = 0, {
                         met <- TCGAanalyze_DMR(data = se,
                                                groupCol = isolate({input$dmrgroupCol}),
                                                group1 = isolate({input$dmrgroup1}),
                                                group2 = isolate({input$dmrgroup2}),
                                                p.cut = isolate({input$dmrpvalue}),
                                                diffmean.cut = isolate({input$dmrthrsld}),
                                                cores = isolate({input$dmrcores}))
                     })
    })
    shinyFileChoose(input, 'dmrfile', roots=volumes, session=session, restrictions=system.file(package='base'))

    dmrdata <- function(){
        inFile <- input$dmrfile
        if (is.null(inFile)) return(NULL)
        file  <- as.character(parseFilePaths(volumes, input$dmrfile)$datapath)
        se <- get(load(file))

        if(class(se)!= class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
            createAlert(session, "dmrmessage", "dmrAlert", title = "Data input error", style =  "danger",
                        content = paste0("Sorry, but I'm expecting a Summarized Experiment object, but I got a: ",
                                         class(se)), append = FALSE)
            return(NULL)
        }
        return(se)

    }

    observeEvent(input$dmrgroupCol , {
        updateSelectizeInput(session, 'dmrgroup1', choices = {
            if (class(dmrdata()) == class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
                if (!is.null(dmrdata()) & input$dmrgroupCol != "" )
                    as.character(colData(dmrdata())[,input$dmrgroupCol])
            }}, server = TRUE)
    })
    observeEvent(input$dmrgroupCol , {
        updateSelectizeInput(session, 'dmrgroup2', choices = {
            if (class(dmrdata()) == class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){

                if (!is.null(dmrdata()) & input$dmrgroupCol != "")
                    as.character(colData(dmrdata())[,input$dmrgroupCol])
            }}, server = TRUE)
    })
    observe({
        updateSelectizeInput(session, 'dmrgroupCol', choices = {
            if(!is.null(dmrdata())) as.character(colnames(colData(dmrdata())))
        }, server = TRUE)
    })
    observe({
        updateSelectizeInput(session, 'meanmetsubgroupCol', choices = {
            if(!is.null(dmrdata())) as.character(colnames(colData(dmrdata())))
        }, server = TRUE)
    })
    observe({
        updateSelectizeInput(session, 'meanmetgroupCol', choices = {
            if(!is.null(dmrdata())) as.character(colnames(colData(dmrdata())))
        }, server = TRUE)
    })


    observeEvent(input$volcanoPlot , {
        output$volcano.plot <- renderPlot({
            if(isolate({input$dmrgroup1}) == "") {
                group1 <- NULL
            } else {
                group1 <- isolate({input$dmrgroup1})
            }

            if(isolate({input$dmrgroup2}) == "") {
                group2 <- NULL
            } else {
                group2 <- isolate({input$dmrgroup2})
            }

            data <- isolate({dmrdata()})
            diffcol <- paste("diffmean",group1,group2,sep = ".")
            pcol <- paste("p.value.adj",group1,group2,sep = ".")


            label <- c("Not Significant",
                       "Hypermethylated",
                       "Hypomethylated")
            label[2:3] <-  paste(label[2:3], "in", group2)
            withProgress(message = 'Creating plot',
                         detail = 'This may take a while...', value = 0, {
                             TCGAVisualize_volcano(x = values(data)[,diffcol],
                                                   y = values(data)[,pcol],
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
                                                   names = NULL,
                                                   x.cut = isolate({as.numeric(input$dmrthrsld)}),
                                                   y.cut = isolate({as.numeric(input$dmrpvalue)}),
                                                   filename = NULL)
                         })

        })})

    observeEvent(input$meanmetPlot , {
        output$mean.plotting <- renderPlot({
            if(isolate({input$meanmetgroupCol}) =="") {
                group <- NULL
            } else {
                group <- isolate({input$meanmetgroupCol})
            }

            if(isolate({input$meanmetsubgroupCol}) =="") {
                subgroup <- NULL
            } else {
                subgroup <- isolate({input$meanmetsubgroupCol})
            }
            withProgress(message = 'Creating plot',
                         detail = 'This may take a while...', value = 0, {
                             TCGAvisualize_meanMethylation(data=dmrdata(),
                                                           groupCol=group,
                                                           subgroupCol=subgroup,
                                                           filename = NULL)
                         })
        })})

    observeEvent(input$meanmetPlot , {
        updateCollapse(session, "collapseDmr", open = "DMR plots")
        output$dmrPlot <- renderUI({
            plotOutput("mean.plotting", width = paste0(isolate({input$meanmetwidth}), "%"), height = isolate({input$meanmetheight}))
        })})
    observeEvent(input$volcanoPlot , {
        updateCollapse(session, "collapseDmr", open = "DMR plots")
        output$dmrPlot <- renderUI({
            plotOutput("volcano.plot", width = paste0(isolate({input$meanmetwidth}), "%"), height = isolate({input$meanmetheight}))
        })})

    output$probesSE <- renderDataTable({
        data <- dmrdata()
        if(!is.null(data)) as.data.frame(values(data))
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


    observeEvent(input$heatmapPlot , {
        output$heatmap.plotting <- renderPlot({
            data <- isolate({dmrdata()})
            colmdata <- isolate({input$colmetadataheatmap})
            rowmdata <- isolate({input$rowmetadataheatmap})
            cluster_rows  <- isolate({input$heatmap.clusterrows})
            show_column_names <- isolate({input$heatmap.show.col.names})
            show_row_names <- isolate({input$heatmap.show.row.names})
            cluster_columns  <- isolate({input$heatmap.clustercol})
            # Get hypo methylated and hypermethylated probes
            idx <- grep("status",colnames(values(data)))
            probes <- which(values(data)[,idx[1]] %in% c("Hypermethylated","Hypomethylated"))
            data <- data[probes,]

            # col.metadata
            col.metadata <- NULL
            print(colmdata)
            if(!is.null(colmdata)) {
                if(length(colmdata) > 0) col.metadata <- subset(colData(data), select=c("patient",colmdata))
            }
            # row.metadata
            row.metadata <- NULL
            print(rowmdata)
            if(!is.null(rowmdata)) {
                if(length(colmdata) > 0) row.metadata <- subset(values(data), select=c(rowmdata))
            }

            withProgress(message = 'Creating plot',
                         detail = 'This may take a while...', value = 0, {
                             p <-  TCGAvisualize_Heatmap(data=assay(data),
                                                         col.metadata=col.metadata,
                                                         row.metadata=row.metadata,
                                                         title = "Heatmap",
                                                         cluster_rows = cluster_rows,
                                                         show_column_names = show_column_names,
                                                         cluster_columns = cluster_columns,
                                                         show_row_names = show_row_names,
                                                         type = "methylation")
                             incProgress(1/2)
                             ComplexHeatmap::draw(p)
                         })
        })})

    observeEvent(input$heatmapPlot , {
        updateCollapse(session, "collapseDmr", open = "DMR plots")
        output$dmrPlot <- renderUI({
            plotOutput("heatmap.plotting", width = paste0(isolate({input$meanmetwidth}), "%"), height = isolate({input$meanmetheight}))
        })})

    observe({
        data <- dmrdata()
        updateSelectizeInput(session, 'colmetadataheatmap', choices = {
            if(!is.null(data)) as.character(colnames(colData(data)))
        }, server = TRUE)
    })
    observe({
        data <- dmrdata()
        updateSelectizeInput(session, 'rowmetadataheatmap', choices = {
            if(!is.null(data)) as.character(colnames(values(data)))
        }, server = TRUE)
    })

    #-------------------- EA
    observeEvent(input$eaplot , {
        updateCollapse(session, "collapseEA", open = "EA plots")
        output$eaPlot <- renderUI({
            plotOutput("ea.plotting", width = paste0(isolate({input$eawidth}), "%"), height = isolate({input$eaheight}))
        })})

    observeEvent(input$eaplot , {
        output$ea.plotting <- renderPlot({
            withProgress(message = 'Creating plot',
                         detail = 'This may take a while...', value = 0, {
                             ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",isolate({input$eagenes}))

                             ResMF <- NULL
                             ResBP <- NULL
                             ResCC <- NULL
                             ResPat <- NULL
                             if(length(grep("NA",ansEA$ResBP)) != ncol(ansEA$ResBP)) ResBP <- ansEA$ResBP
                             if(length(grep("NA",ansEA$ResCC)) != ncol(ansEA$ResCC)) ResCC <- ansEA$ResCC
                             if(length(grep("NA",ansEA$ResMF)) != ncol(ansEA$ResMF)) ResMF <- ansEA$ResMF
                             if(length(grep("NA",ansEA$ResPat)) != ncol(ansEA$ResPat)) ResPat <- ansEA$ResPat

                             # Enrichment Analysis EA (TCGAVisualize)
                             # Gene Ontology (GO) and Pathway enrichment barPlot

                             TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),
                                                     GOBPTab = ResBP,
                                                     GOCCTab = ResCC,
                                                     GOMFTab = ResMF,
                                                     PathTab = ResPat,
                                                     color = c(isolate({input$colBP}),
                                                               isolate({input$colCC}),
                                                               isolate({input$colMF}),
                                                               isolate({input$colPat})),
                                                     nRGTab = isolate({input$eagenes}),
                                                     nBar = isolate({input$nBar}),
                                                     filename = NULL)
                         })
        })})

}
