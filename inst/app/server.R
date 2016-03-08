library(shiny)
library(shinyFiles)
library(SummarizedExperiment)
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
    output$directorypath <- renderText({parseDirPath(volumes, input$folder)})
    #----------------- Ontology
    output$ontSearchLink <-  renderText({
        getPath <- parseDirPath(volumes, input$folder)
        getFtype <- input$ontftypeFilter
        if (input$ontSearchDownloadBt ) {
            link <- c()

            df <- data.frame(database = input$allRows[seq(1, length(input$allRows), 5)],
                             ID = input$allRows[seq(2, length(input$allRows), 5)],
                             Sample = input$allRows[seq(3, length(input$allRows), 5)],
                             Experiment = input$allRows[seq(4, length(input$allRows), 5)],
                             organism = input$allRows[seq(5, length(input$allRows), 5)])
            if(length(getPath) > 0) {
                biOmicsDownload(df, path = getPath, enc.file.type = getFtype, rmap.file.type = getFtype)
            } else {
                biOmicsDownload(df, enc.file.type = getFtype, rmap.file.type = getFtype)
            }
            print("End of download")
        }}
    )

    dataInput <- reactive({
        indexes <- c()
        query <- data.frame()
        link <- c()
        term <- input$ontSamplesFilter
        exp <-  input$ontExpFilter
        query <- biOmicsSearch(input$ontSamplesFilter,
                               experiment = input$ontExpFilter)
        # improve using subset - subset(data,selection,projection)
        biosample.encode <- get.obj("biosample.encode")
        return(list(system = names(sort(table(biosample.encode[biosample.encode$biosample %in% query$Sample,]$system),decreasing = T)[1]),
                    result = query))
    })
    observe({
        if(input$ontReport){
            create.report(dataInput()$result, system = dataInput()$system)
            create.report(dataInput()$result, path = paste0(getwd(),"/report") ,system = dataInput()$system)
        }
    })
    output$system <-  renderText({
        if (input$ontSearchBt) {
            paste0("System: ", dataInput()$system)
        }
    })

    output$tcgaprepare <-  renderText({
        if (input$tcgaPrepareBt) {
            # read the data from the downloaded path
            # prepare it

            # Files types
            ftype <- NULL
            rnaseqFtype <- input$tcgaFrnaseqtypeFilter
            rnaseqv2Ftype <- input$tcgaFrnaseqv2typeFilter
            gwsFtype <- input$tcgaFgwstypeFilter

            # Dir to saved the files
            getPath <- parseDirPath(volumes, input$tcgafolder)
            if (length(getPath) == 0) getPath <- "."
            samplesType <- input$tcgasamplestypeFilter

            filename <- parseDirPath(volumes, input$tcgapreparefolder)
            if(length(filename) == 0) filename <- "."
            if(length(input$tcgapreparefile) > 0) filename <- file.path(filename,input$tcgapreparefile)

            withProgress(message = 'Prepare in progress',
                         detail = 'This may take a while...', value = 0, {
                             df <- data.frame(name = input$allRows[seq(6, length(input$allRows), 7)])
                             x <- TCGAquery()
                             x <- x[x$name %in% df$name,]
                             for (i in unique(x$Platform)) {
                                 incProgress(1/length(unique(x$Platform)))
                                 aux <- x[x$Platform ==i,]
                                 ftype <- NULL
                                 if (i == "IlluminaHiSeq_RNASeqV2") ftype <- rnaseqv2Ftype
                                 if (i == "IlluminaHiSeq_RNASeq") ftype <- rnaseqFtype
                                 if (i == "Genome_Wide_SNP_6") ftype <- gwsFtype
                                 if (length(samplesType) == 0) {
                                     samples <- NULL
                                 } else {
                                     samples <- unlist(lapply(samplesType,function(type){
                                         s <- unlist(str_split(aux$barcode,","))
                                         s[grep(type,substr(s,14,15))]
                                     }))
                                 }
                                 if(!(length(unique(x$Platform)) == 1 &
                                      length(input$tcgapreparefile) > 0)){
                                     filename <- NULL
                                 }
                                 trash <- TCGAprepare(aux,dir = getPath,
                                                      summarizedExperiment = as.logical(input$prepareRb),
                                                      save = TRUE,
                                                      type = ftype,
                                                      filename=filename,
                                                      samples = samples)
                             }

                         })
            print("End of Prepare")
        }
    })


    output$ontSearchtbl <- renderDataTable({
        if (input$ontSearchBt) {
            dataInput()$result
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
    )

    #------------- TCGA ------------------

    # Table render
    output$tcgaSearchtbl <- renderDataTable({
        tumor = input$tcgaTumorFilter
        platform = input$tcgaExpFilter
        level = input$tcgaLevelFilter
        if (input$tcgaSearchBt) {
            tbl <- TCGAquery(tumor = tumor,
                             platform = platform,
                             level = level)
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
    )

    # Download
    shinyDirChoose(input, 'tcgafolder', roots=volumes, session=session,
                   restrictions=system.file(package='base'))
    output$tcgadirectorypath <- renderText({parseDirPath(volumes, input$tcgafolder)})
    shinyDirChoose(input, 'tcgapreparefolder', roots=volumes, session=session,
                   restrictions=system.file(package='base'))
    output$tcgapreparedir <- renderText({parseDirPath(volumes, input$tcgapreparefolder)})


    output$tcgaSearchLink <-  renderText({

        # Files types
        ftype <- NULL
        rnaseqFtype <- input$tcgaFrnaseqtypeFilter
        rnaseqv2Ftype <- input$tcgaFrnaseqv2typeFilter
        gwsFtype <- input$tcgaFgwstypeFilter

        # Dir to save the files
        getPath <- parseDirPath(volumes, input$tcgafolder)
        if (length(getPath) == 0) getPath <- "."
        samplesType <- input$tcgasamplestypeFilter
        if (input$tcgaDownloadBt ) {

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
            print("End of download")
        }

    })

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
    output$mafdirectorypath <- renderText({parseDirPath(volumes, input$maffolder)})

    output$mafDownloadfTxt <-  renderText({
        maf.files <- get.obj("maf.files")
        # Dir to save the files
        getPath <- parseDirPath(volumes, input$maffolder)
        if (length(getPath) == 0) getPath <- "."
        if (input$mafDownloadBt ) {

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
            print("End of download")
        }

    })
    mut <- function(){
        inFile <- input$maffile

        if (is.null(inFile)) return(NULL)
        ret <- read.table(inFile$datapath, fill = TRUE,
                          comment.char = "#", header = TRUE, sep = "\t", quote='')
        return(ret)

    }
    plotting <- reactive({
        if(input$oncoprintPlot){
            mut <- mut()
            create.oncoprint(mut=mut,genes=input$oncoGenes)
        }
    })
    output$oncoPlot <- renderPlot({
        if(input$oncoprintPlot){
            plotting()
        }
    })
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
    output$tcgadmr <-  renderText({
        if (input$dmrAnalysis) {
            # read the data from the downloaded path
            # prepare it
            se <- dmrdata()
            se <- subset(se,subset = (rowSums(is.na(assay(se))) == 0))
            withProgress(message = 'DME analysis in progress',
                         detail = 'This may take a while...', value = 0, {
                             met <- TCGAanalyze_DMR(data = se,
                                             groupCol = input$dmrgroupCol,
                                             group1 = input$dmrgroup1,
                                             group2 = input$dmrgroup2,
                                             p.cut=input$dmrpvalue,
                                             diffmean.cut=input$dmrthrsld,
                                             cores = input$dmrcores)

                         })
            print("End of DMR analysis")
        }
    })

    dmrdata <- function(){
        inFile <- input$dmrfile
        if (is.null(inFile)) return(NULL)
        se <- get(load(input$dmrfile$datapath))
        return(se)
    }

    observe({
        updateSelectizeInput(session, 'dmrgroup1', choices = {
            if(!is.null(dmrdata()) & input$dmrgroupCol !="")
                as.character(colData(dmrdata())[,input$dmrgroupCol])
        }, server = TRUE)
    })
    observe({
        updateSelectizeInput(session, 'dmrgroup2', choices = {

            if(!is.null(dmrdata()) & input$dmrgroupCol !="")
                as.character(colData(dmrdata())[,input$dmrgroupCol])
            }, server = TRUE)
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

}

