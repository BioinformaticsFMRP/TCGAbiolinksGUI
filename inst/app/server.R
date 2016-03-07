library(shiny)
library(shinyFiles)
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
        return(list(system = names(sort(table(biosample.encode[biosample.encode$biosample %in% query$Sample,]$system),decreasing = T)[1]),
                    result = query))
    })

    output$system <-  renderText({
        if (input$ontSearchBt) {
            paste0("System: ", dataInput()$system)
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
                                 if (x$Platform == "IlluminaHiSeq_RNASeqV2") ftype <- rnaseqv2Ftype
                                 if (x$Platform == "IlluminaHiSeq_RNASeq") ftype <- rnaseqFtype
                                 if (x$Platform == "Genome_Wide_SNP_6") ftype <- gwsFtype
                                 if (length(samplesType) == 0) {
                                     samples <- NULL
                                 } else {
                                     samples <- unlist(lapply(samplesType,function(type){
                                         s <- unlist(str_split(x$barcode,","))
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
        all.df[,c(10,1,2,4,5,7)]
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

}

