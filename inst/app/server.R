library(shiny)
library(shinyFiles)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(shinyBS)
library(stringr)
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
                                 print(samples)
                             }
                             if(is.null(samples)){
                                 TCGAdownload(x, path = getPath,type = ftype)
                             } else if (length(samples) > 0){
                                 TCGAdownload(x, path = getPath,type = ftype,samples = samples)
                             }
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
    #----------------------------------------------------------
    #  Profile plot
    #-----------------------------------------------------------
    observe({
        data <- profileplotdata()
        updateSelectizeInput(session, 'profileplotgroup', choices = {
            if(!is.null(data)) as.character(colnames(data))
        }, server = TRUE)
    })

    observe({
        data <- profileplotdata()
        updateSelectizeInput(session, 'profileplotsubtype', choices = {
            if(!is.null(data)) as.character(colnames(data))
        }, server = TRUE)
    })

    profileplotdata <- function(){
        inFile <- input$profileplotfile
        if (is.null(inFile)) return(NULL)
        file  <- as.character(parseFilePaths(volumes, inFile)$datapath)
        df <- get(load(file))

        if (class(df) != class(data.frame())){
            createAlert(session, "dmrmessage", "dmrAlert", title = "Data input error", style =  "danger",
                        content = paste0("Sorry, but I'm expecting a data frameobject, but I got a: ",
                                         class(df)), append = FALSE)
            return(NULL)
        }
        return(df)
    }

    observe({
        groupCol <-  input$profileplotgroup
        na.rm.groups <-  input$profileplotrmnagroup
        data <- isolate({profileplotdata()})


        if(na.rm.groups){
            data <- data[!is.na(data[,groupCol]),]
            data <- data[which(data[,groupCol] != "NA"),]
        }
        x <- length(unique(data[,groupCol]))
        m1 <- 0
        m3 <- 0

        if (x == 2) { m1 <- -5.0; m3 <-  1.8}
        if (x == 3) { m1 <- -5.0; m3 <-  0.0}
        if (x == 5) { m1 <- -4.8; m3 <- -1.5}
        if (x == 6) { m1 <- -4.8; m3 <- -2.0}
        if (x == 7) { m1 <- -4.8; m3 <- -2.1}
        if (x == 8) { m1 <- -4.8; m3 <- -2.5}

        # Control the value, min, max, and step.
        # Step size is 2 when input value is even; 1 when value is odd.
        updateSliderInput(session, "margin1", value = m1,
                          min = -10, max = 10, step = 0.1)
        updateSliderInput(session, "margin3", value = m3,
                          min = -10, max = 10, step = 0.1)

    })

    observeEvent(input$profileplotBt , {
        output$profile.plotting <- renderPlot({
            data <- isolate({profileplotdata()})
            subtypeCol <- isolate({input$profileplotsubtype})
            groupCol <-  isolate({input$profileplotgroup})
            na.rm.groups <-  isolate({input$profileplotrmnagroup})
            na.rm.subtypes  <-  isolate({input$profileplotrmnasub})
            m1  <-  isolate({input$margin1})
            m2  <-  isolate({input$margin2})
            m3  <-  isolate({input$margin3})
            m4  <-  isolate({input$margin4})

            withProgress(message = 'Creating plot',
                         detail = 'This may take a while...', value = 0, {

                             TCGAvisualize_profilePlot(data = data,
                                                       groupCol=groupCol,
                                                       subtypeCol=subtypeCol,
                                                       na.rm.groups = na.rm.groups,
                                                       na.rm.subtypes = na.rm.subtypes,
                                                       plot.margin=c(m1,m2,m3,m4))

                         })
        })})

    observeEvent(input$profileplotBt , {
        updateCollapse(session, "collapseprofileplot", open = "Profile plot")
        output$profileplot <- renderUI({
            plotOutput("profile.plotting", width = paste0(isolate({input$profilewidth}), "%"), height = isolate({input$profileheight}))
        })})
    shinyFileChoose(input, 'profileplotfile', roots=volumes, session=session, restrictions=system.file(package='base'))

    #------------------------------------------------
    # Survival plot
    # -----------------------------------------------
    #----------------------------------------------------------
    #  Profile plot
    #-----------------------------------------------------------
    observe({
        data <- survivalplotdata()
        updateSelectizeInput(session, 'survivalplotgroup', choices = {
            if(!is.null(data)) as.character(colnames(data))
        }, server = TRUE)
    })

    observe({
        data <- survivalplotdata()
        updateSelectizeInput(session, 'survivalplotsubtype', choices = {
            if(!is.null(data)) as.character(colnames(data))
        }, server = TRUE)
    })

    survivalplotdata <- function(){
        inFile <- input$survivalplotfile
        if (is.null(inFile)) return(NULL)
        file  <- as.character(parseFilePaths(volumes, inFile)$datapath)
        df <- get(load(file))

        if (class(df) != class(data.frame())){
            closeAlert(session, "survivalAlert")
            createAlert(session, "survivalmessage", "survivalAlert", title = "Data input error", style =  "danger",
                        content = paste0("Sorry, but I'm expecting a data frame object, but I got a: ",
                                         class(df)), append = FALSE)
            return(NULL)
        }
        return(df)
    }

    observeEvent(input$survivalplotBt , {
        output$survival.plotting <- renderPlot({
            data <- isolate({survivalplotdata()})
            clusterCol <-  isolate({input$survivalplotgroup})
            closeAlert(session, "survivalAlert")
            if(length(unique(data[,clusterCol])) == 1){
                createAlert(session, "survivalmessage", "survivalAlert", title = "Data input error", style =  "danger",
                            content = paste0("Sorry, but I'm expecting at least two groups"), append = FALSE)
                return(NULL)
            }
            withProgress(message = 'Creating plot',
                         detail = 'This may take a while...', value = 0, {

                             TCGAanalyze_survival(data = data,
                                                  clusterCol = clusterCol,
                                                  filename = NULL)

                         })
        })})

    observeEvent(input$survivalplotBt , {
        updateCollapse(session, "collapsesurvivalplot", open = "survival plot")
        output$survivalplot <- renderUI({
            plotOutput("survival.plotting", width = paste0(isolate({input$survivalwidth}), "%"), height = isolate({input$survivalheight}))
        })})
    shinyFileChoose(input, 'survivalplotfile', roots=volumes, session=session, restrictions=system.file(package='base'))

    # -------------------------------------------------
    # DEA
    # -------------------------------------------------
    observeEvent(input$deaAnalysis , {
        # read the data from the downloaded path
        # prepare it
        se <- isolate({deadata()})

        g1 <- isolate({input$deagroup1})
        g2 <- isolate({input$deagroup2})
        groupCol <-  isolate({input$deagroupCol})
        idx.g1 <- which(colData(se)[,groupCol] == g1)
        samples.g1 <- colData(se)[idx.g1,"barcode"]
        idx.g2 <- which(colData(se)[,groupCol] == g2)
        samples.g2 <- colData(se)[idx.g2,"barcode"]
        method <- isolate({input$deamethod})
        fdr.cut <- isolate({input$deapvalue})
        logFC.cut <- isolate({input$deathrsld})
        withProgress(message = 'dea analysis in progress',
                     detail = 'This may take a while...', value = 0, {


                         # normalization of genes
                         if(isolate({input$deanormalization})) {
                             exp <- TCGAanalyze_Normalization(tabDF = assay(se),
                                                              geneInfo = TCGAbiolinks::geneInfo,
                                                              method = isolate({input$deanormalizationmet})
                             )
                         }
                         # quantile filter of genes
                         if(isolate({input$deafilter})) {
                             dataFilt <- TCGAanalyze_Filtering(tabDF = exp,
                                                               method = isolate({input$deafilteringmet}),
                                                               qnt.cut =  isolate({input$deafilteringcut}))
                         }

                         exp <- TCGAanalyze_DEA(mat1 = assay(se[,samples.g1]),
                                                mat2 = assay(se[,samples.g2]),
                                                Cond1type = g1 ,
                                                Cond2type = g2,
                                                #fdr.cut  = fdr.cut,
                                                #logFC.cut = logFC.cut,
                                                method = method)
                         exp <- TCGAanalyze_LevelTab(exp,
                                                     typeCond1 = g1,
                                                     typeCond2 = g2,
                                                     TableCond1 = assay(se[,samples.g1]),
                                                     TableCond2 = assay(se[,samples.g2]))
                         exp$status <- "Insignificant"
                         exp[exp$logFC >= logFC.cut & exp$FDR <= fdr.cut,"status"] <- paste0("Upregulated in ", g1)
                         exp[exp$logFC <= -logFC.cut & exp$FDR <= fdr.cut,"status"] <- paste0("Downregulated in ", g1)
                     })


        save(exp, file = paste0("result_dea_", g1, "_", g2,".rda"))
        createAlert(session, "deamessage", "deaAlert", title = "DEA completed", style =  "danger",
                    content = paste0("Results saved in: result_dea_", g1, "_", g2,".rda"), append = FALSE)
    })
    shinyFileChoose(input, 'deafile', roots=volumes, session=session, restrictions=system.file(package='base'))

    deadata <- function(){
        inFile <- input$deafile
        if (is.null(inFile)) return(NULL)
        file  <- as.character(parseFilePaths(volumes, input$deafile)$datapath)
        se <- get(load(file))

        if(class(se)!= class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
            createAlert(session, "deamessage", "deaAlert", title = "Data input error", style =  "danger",
                        content = paste0("Sorry, but I'm expecting a Summarized Experiment object, but I got a: ",
                                         class(se)), append = FALSE)
            return(NULL)
        }
        return(se)
    }

    observeEvent(input$deagroupCol , {
        updateSelectizeInput(session, 'deagroup1', choices = {
            if (class(deadata()) == class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
                if (!is.null(deadata()) & input$deagroupCol != "" )
                    as.character(colData(deadata())[,input$deagroupCol])
            }}, server = TRUE)
    })
    observeEvent(input$deagroupCol , {
        updateSelectizeInput(session, 'deagroup2', choices = {
            if (class(deadata()) == class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){

                if (!is.null(deadata()) & input$deagroupCol != "")
                    as.character(colData(deadata())[,input$deagroupCol])
            }}, server = TRUE)
    })
    observe({
        updateSelectizeInput(session, 'deagroupCol', choices = {
            if(!is.null(deadata())) as.character(colnames(colData(deadata())))
        }, server = TRUE)
    })

    observeEvent(input$volcanodeaPlot , {
        output$volcano.dea.plot <- renderPlot({
            g1 <- isolate({input$deagroup1})
            g2 <- isolate({input$deagroup2})
            fdr.cut <- isolate({input$deapvalue})
            logFC.cut <- isolate({input$deathrsld})
            dea.result <- get(load( paste0("result_dea_", g1, "_", g2,".rda")))

            label <- c("Not Significant",
                       "Upregulated",
                       "Downregulated")
            label[2:3] <-  paste(label[2:3], "in", g2)


            withProgress(message = 'Creating plot',
                         detail = 'This may take a while...', value = 0, {
                             TCGAVisualize_volcano(x = dea.result$logFC,
                                                   y = dea.result$FDR,
                                                   ylab =   expression(paste(-Log[10],
                                                                             " (FDR corrected -P values)")),
                                                   xlab = " Gene expression fold change (Log2)",
                                                   color = c(isolate({input$coldeainsignificant}),
                                                             isolate({input$colUpregulated}),
                                                             isolate({input$colDownregulated})),
                                                   title =  paste("Volcano plot", "(", g2, "vs", g1,")"),
                                                   legend=  "Legend",
                                                   label = label,
                                                   names = NULL,
                                                   x.cut = isolate({as.numeric(input$deathrsld)}),
                                                   y.cut = isolate({as.numeric(input$deapvalue)}),
                                                   filename = NULL)
                         })

        })})

    observeEvent(input$volcanodeaPlot , {
        updateCollapse(session, "collapsedea", open = "dea plots")
        output$deaPlot <- renderUI({
            plotOutput("volcano.dea.plot", width = paste0(isolate({input$deawidth}), "%"), height = isolate({input$deaheight}))
        })})

    output$deaSE <- renderDataTable({
        data <- deadata()
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

    #----------------------------------------------
    #                 DEA Pathview
    observeEvent(input$pathwaygraphBt , {

        pathway.id <- isolate({input$pathway.id})
        kegg.native <- isolate({input$kegg.native.checkbt})
        g1 <- isolate({input$deagroup1})
        g2 <- isolate({input$deagroup2})
        data <- get(load( paste0("result_dea_", g1, "_", g2,".rda")))

        gene <- strsplit(data$mRNA,"\\|")
        data$SYMBOL <- unlist(lapply(gene,function(x) x[1]))

        # Converting Gene symbol to geneID
        library(clusterProfiler)
        eg = as.data.frame(bitr(data$SYMBOL,
                                fromType="SYMBOL",
                                toType="ENTREZID",
                                annoDb="org.Hs.eg.db"))
        eg <- eg[!duplicated(eg$SYMBOL),]

        data <- merge(data,eg,by="SYMBOL")

        data <- subset(data, select = c("ENTREZID", "logFC"))
        genelistDEGs <- as.numeric(data$logFC)
        names(genelistDEGs) <- data$ENTREZID
        withProgress(message = 'Creating pathway graph',
                     detail = 'This may take a while...', value = 0, {
                         require("pathview")
                         # pathway.id: hsa05214 is the glioma pathway
                         # limit: sets the limit for gene expression legend and color
                         hsa05214 <- pathview(gene.data  = genelistDEGs,
                                              pathway.id = pathway.id,
                                              species    = "hsa",
                                              kegg.native = kegg.native,
                                              limit      = list(gene=as.integer(max(abs(genelistDEGs)))))
                     })
        if(kegg.native) {
            extension <- ".pathview.png"
        } else {
            extension <- ".pathview.pdf"
        }

        createAlert(session, "deamessage", "deaAlert", title = "Pathway graph created", style =  "active",
                    content = paste0("Results saved in: ", pathway.id,extension), append = FALSE)

    })
    #----------------------------------------------
    starburst <- function(){
        g1 <- isolate({input$starburstgroup1})
        g2 <- isolate({input$starburstgroup2})
        logFC.cut <- isolate({input$starburstexpFC})
        exp.p.cut <- isolate({input$starburstexFDR})
        diffmean.cut <- isolate({input$starburstmetdiff})
        met.p.cut <- isolate({input$starburstmetFDR})
        exp <- result.dea.data()
        met <-result.dmr.data()

        colors <- c(isolate({input$sbcolInsignigicant}),
                    isolate({input$sbcolUpHypo}),
                    isolate({input$sbcolDownHypo}),
                    isolate({input$sbcolHypo}),
                    isolate({input$sbcolHyper}),
                    isolate({input$sbcolUp}),
                    isolate({input$sbcolDown}),
                    isolate({input$sbcolUpHyper}),
                    isolate({input$sbcolDownHyper}))
        result <- TCGAvisualize_starburst(met = met,
                                          exp = exp,
                                          group1 = g1,
                                          group2 = g2,
                                          color = colors,
                                          exp.p.cut = exp.p.cut,
                                          met.p.cut = met.p.cut,
                                          diffmean.cut = diffmean.cut,
                                          logFC.cut = logFC.cut,
                                          return.plot = TRUE)
    }
    # -------------- Starburst plot
    observeEvent(input$starburstPlot , {
        output$starburst.plot <- renderPlot({
            withProgress(message = 'Creating plot',
                         detail = 'This may take a while...', value = 0, {
                             starburst()$plot
                         })
        })})

    observeEvent(input$starburstPlot , {
        updateCollapse(session, "collapsedea", open = "dea plots")
        output$starburstPlot <- renderUI({
            plotOutput("starburst.plot", width = paste0(isolate({input$starburstwidth}), "%"), height = isolate({input$starburstheight}))
        })})
    observe({
        updateSelectizeInput(session, 'starburstgroup1', choices = {
            if(!is.null(result.dea.data())) {
                x <- as.character(colnames(result.dea.data()))
                x[-which(x %in% c("mRNA", "logFC","FDR", "Delta","status"))]
            }
        }, server = TRUE)
    })
    observe({
        updateSelectizeInput(session, 'starburstgroup2', choices = {
            if(!is.null(result.dea.data())) {
                x <- as.character(colnames(result.dea.data()))
                x[-which(x %in% c("mRNA", "logFC","FDR", "Delta","status"))]
            }
        }, server = TRUE)
    })
    result.dea.data <- function(){
        inFile <- input$starburstexpfile
        if (is.null(inFile)) return(NULL)
        file  <- as.character(parseFilePaths(volumes, input$starburstexpfile)$datapath)
        se <- get(load(file))

        if(class(se)!= class(data.frame())){
            createAlert(session, "deamessage", "deaAlert", title = "Data input error", style =  "danger",
                        content = paste0("Sorry, but I'm expecting a Data frame object, but I got a: ",
                                         class(se)), append = FALSE)
            return(NULL)
        }
        return(se)
    }
    result.dmr.data <- function(){
        inFile <- input$starburstmetfile
        if (is.null(inFile)) return(NULL)
        file  <- as.character(parseFilePaths(volumes, input$starburstmetfile)$datapath)
        se <- get(load(file))

        if(class(se)!= class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
            createAlert(session, "deamessage", "deaAlert", title = "Data input error", style =  "danger",
                        content = paste0("Sorry, but I'm expecting a Summarized Experiment object, but I got a: ",
                                         class(se)), append = FALSE)
            return(NULL)
        }
        return(se)
    }
    shinyFileChoose(input, 'starburstmetfile', roots=volumes, session=session, restrictions=system.file(package='base'))
    shinyFileChoose(input, 'starburstexpfile', roots=volumes, session=session, restrictions=system.file(package='base'))


    output$starburstResult <- renderDataTable({
        data <- starburst()$starburst
        if(!is.null(data)) as.data.frame(data)
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

}
