suppressPackageStartupMessages({
    library(shiny)
    library(shinyFiles)
    library(SummarizedExperiment)
    library(MultiAssayExperiment)
    library(TCGAbiolinks)
    library(ggplot2)
    library(shinyBS)
    library(stringr)
    library(ggrepel)
    library(plotly)
    library(pathview)
    library(htmlwidgets)
    library(ELMER)
    library(readr)
    library(data.table)
    library(grid)
    library(maftools)
    library(dplyr)
    data(maf.tumor)
    data(GDCdisease)
    options(shiny.maxRequestSize=-1) # Remove limit of upload
    options(shiny.deprecation.messages=FALSE)
    options(warn =-1)
})

getDataCategory <- function(legacy){
    data.category.hamonirzed <- sort(c("Transcriptome Profiling",
                                       "Copy Number Variation",
                                       "Simple Nucleotide Variation",
                                       "DNA Methylation"))

    data.category.legacy <- sort(c("Copy number variation",
                                   "Simple Nucleotide Variation",
                                   "Raw Sequencing Data",
                                   "Protein expression",
                                   "Gene expression",
                                   "DNA methylation",
                                   "Raw microarray data",
                                   # "Structural Rearrangement", # Controlled
                                   "Other"))
    if(legacy) return(data.category.legacy)
    return(data.category.hamonirzed)
}

createTable <- function(df,tableType = "TCGAbiolinks"){
    DT::datatable(df,
                  extensions = c('Buttons',"FixedHeader"),
                  class = 'cell-border stripe',
                  options = list(dom = 'Blfrtip',
                                 buttons =
                                     list('colvis', list(
                                         extend = 'collection',
                                         buttons = list(list(extend='csv',
                                                             filename = tableType),
                                                        list(extend='excel',
                                                             filename = tableType),
                                                        list(extend='pdf',
                                                             title = "",
                                                             filename= tableType)),
                                         text = 'Download'
                                     )),
                                 fixedHeader = TRUE,
                                 pageLength = 20,
                                 scrollX = TRUE,
                                 lengthMenu = list(c(10, 20, -1), c('10', '20', 'All'))
                  ),
                  filter = 'top'
    )
}

getFileType <-  function(legacy, data.category){
    file.type <- NULL
    if(grepl("Copy number variation",data.category, ignore.case = TRUE) & legacy)
        file.type <- c("nocnv_hg18.seg",
                       "nocnv_hg19.seg",
                       "hg19.seg",
                       "hg18.seg")
    if(grepl("Gene expression", data.category, ignore.case = TRUE)  & legacy)
        file.type <- c("normalized_results",
                       "results")

    if(data.category == "Raw microarray data") file.type <- c("idat")

    return(file.type)
}

getExpStrategy <-  function(legacy, platform){
    experimental.strategy <- NULL

    # These are the cases we need to distinguish
    if(grepl("Illumina HiSeq", platform, ignore.case = TRUE)  & legacy)
        experimental.strategy <- c("Total RNA-Seq",
                                   "RNA-Seq",
                                   "miRNA-Seq")
    if(grepl("Illumina GA", platform, ignore.case = TRUE)  & legacy)
        experimental.strategy <- c("RNA-Seq",
                                   "miRNA-Seq")

    return(experimental.strategy)
}


getWorkFlow <-  function(legacy, data.category){
    workflow <- NULL
    if(data.category == "Transcriptome Profiling" & !legacy)
        workflow <- c("HTSeq - Counts",
                      "HTSeq - FPKM-UQ",
                      "HTSeq - FPKM")
    return(workflow)
}

getPlatform <-  function(legacy, data.category){
    platform <- NULL

    if(!legacy & data.category != "DNA Methylation" ) return(platform) # platform is not used for harmonized
    if(grepl("Copy number variation",data.category, ignore.case = TRUE)) platform <- "Affymetrix SNP Array 6.0"
    if(data.category == "Protein expression") platform <- "MDA RPPA Core"
    if(data.category == "Gene expression") platform <- c("Illumina HiSeq",
                                                         "HT_HG-U133A",
                                                         "AgilentG4502A_07_2",
                                                         "AgilentG4502A_07_1",
                                                         "HuEx-1_0-st-v2")
    if(data.category == "DNA methylation") platform <- c("Illumina Human Methylation 450",
                                                         "Illumina Human Methylation 27",
                                                         "Illumina DNA Methylation OMA003 CPI",
                                                         "Illumina DNA Methylation OMA002 CPI",
                                                         "Illumina Hi Seq")
    if(data.category == "DNA Methylation") platform <- c("Illumina Human Methylation 450",
                                                         "Illumina Human Methylation 27")

    if(data.category == "Raw microarray data") platform <- c("Illumina Human Methylation 450",
                                                             "Illumina Human Methylation 27")
    return(platform)
}

getDataType <- function(legacy, data.category){
    data.type <- NULL
    if(data.category == "Transcriptome Profiling" & !legacy)
        data.type <- c("Gene Expression Quantification",
                       "Isoform Expression Quantification",
                       "miRNA Expression Quantification")

    if(grepl("Copy number variation",data.category, ignore.case = TRUE) & !legacy)
        data.type <- c("Copy Number Segment",
                       "Masked Copy Number Segment")

    if(data.category == "Gene expression" & !legacy)
        data.type <- c("Gene Expression Quantification",
                       "Isoform Expression Quantification",
                       "miRNA Expression Quantification")

    if(data.category == "Gene expression" & legacy)
        data.type  <- c("Gene expression quantification",
                        "miRNA gene quantification",
                        "Exon junction quantification",
                        "Exon quantification",
                        "miRNA isoform quantification")
    if(data.category == "Raw microarray data" & legacy)
        data.type  <- c("Raw intensities")

    return(data.type)
}

table.code <- c('01','02','03','04','05','06','07','08','09','10',
                '11','12','13','14','20','40','50','60','61')
names(table.code) <- c("Primary solid Tumor","Recurrent Solid Tumor",
                       "Primary Blood Derived Cancer - Peripheral Blood",
                       "Recurrent Blood Derived Cancer - Bone Marrow",
                       "Additional - New Primary",
                       "Metastatic","Additional Metastatic",
                       "Human Tumor Original Cells",
                       "Primary Blood Derived Cancer - Bone Marrow",
                       "Blood Derived Normal","Solid Tissue Normal",
                       "Buccal Cell Normal","EBV Immortalized Normal",
                       "Bone Marrow Normal","Control Analyte",
                       "Recurrent Blood Derived Cancer - Peripheral Blood",
                       "Cell Lines","Primary Xenograft Tissue",
                       "Cell Line Derived Xenograft Tissue")

tcga.code <- c("Primary solid Tumor","Recurrent Solid Tumor",
               "Primary Blood Derived Cancer - Peripheral Blood",
               "Recurrent Blood Derived Cancer - Bone Marrow",
               "Additional - New Primary",
               "Metastatic","Additional Metastatic",
               "Human Tumor Original Cells",
               "Primary Blood Derived Cancer - Bone Marrow",
               "Blood Derived Normal","Solid Tissue Normal",
               "Buccal Cell Normal","EBV Immortalized Normal",
               "Bone Marrow Normal","Control Analyte",
               "Recurrent Blood Derived Cancer - Peripheral Blood",
               "Cell Lines","Primary Xenograft Tissue",
               "Cell Line Derived Xenograft Tissue")
names(tcga.code) <- c('01','02','03','04','05','06','07','08','09','10',
                      '11','12','13','14','20','40','50','60','61')

getMatchedPlatform <- function(query){
    matched <- NULL
    for(plat in query$Platform){
        aux <- query[query$Platform == plat,]
        if(is.null(matched)){
            matched <- unlist(str_split(aux$barcode,","))
            matched <- substr(matched,1,15)
        } else {
            barcode <- unlist(str_split(aux$barcode,","))
            barcode <- substr(barcode,1,15)
            matched <- intersect(matched, barcode)
        }
    }
    return(matched)
}


getMatchedType <- function(barcode,type){

    code <- c("TP","TR","TB","TRBM","TAP","TM","TAM","THOC",
              "TBM","NB","NT","NBC","NEBV","NBM","CELLC","TRB",
              "CELL","XP","XCL")

    names(code) <- c('01','02','03','04','05','06','07','08','09','10',
                     '11','12','13','14','20','40','50','60','61')

    type <- code[type]
    groups <- t(combn(type,2))
    matched <- NULL
    for(i in 1:nrow(groups)) {
        if(is.null(matched)){
            matched <- TCGAquery_MatchedCoupledSampleTypes(unique(barcode),
                                                           c(groups[i,1], groups[i,2]))
            matched <- substr(matched,1,15)
        } else {
            aux <- TCGAquery_MatchedCoupledSampleTypes(unique(barcode),
                                                       c(groups[i,1], groups[i,2]))
            aux <- substr(aux,1,15)
            matched <- intersect(matched, aux)
        }
    }
    return(matched)
}


# This will be used to parse the text areas input
# possibilities of separation , ; \n
parse.textarea.input <- function(text){
    sep <- NULL
    if(grepl(";",text)) sep <- ";"
    if(grepl(",",text)) sep <- ","
    if(grepl("\n",text)) sep <- "\n"
    if(is.null(sep)) {
        text <- text
    } else {
        text <- unlist(stringr::str_split(text,sep))
    }
    return (text)
}

get.volumes <- function(directory = NULL){

    if(is.null(directory) ||  identical(directory, character(0))) {
        volumes <- c(wd="./",home=Sys.getenv("HOME"), getVolumes()(), temp=tempdir())
    } else {
        volumes <- c(home=Sys.getenv("HOME"), getVolumes()(), temp=tempdir(),wd="./")
        path <- parseDirPath(volumes, directory)
        names(path) <- path
        volumes <- c(path, home=Sys.getenv("HOME"), getVolumes()(), temp=tempdir(),wd="./")
    }
    return(volumes)
}

#' @title  Server side
#' @description Server side
#' @param input - input signal
#' @param output - output signal
#' @importFrom downloader download
#' @import pathview ELMER TCGAbiolinks SummarizedExperiment shiny ggrepel UpSetR
#' @keywords internal
TCGAbiolinksGUIServer <- function(input, output, session) {

    session$onSessionEnded(stopApp)
    server.path <- ifelse(system.file("app", package = "TCGAbiolinksGUI") == "",
                          "server",
                          file.path(system.file("app", package = "TCGAbiolinksGUI"),"server"))

    source(file.path(server.path, "getmolecular.R"),  local = TRUE)$value
    source(file.path(server.path, "getsubtype.R"),  local = TRUE)$value
    source(file.path(server.path, "getmutation.R"),  local = TRUE)$value
    source(file.path(server.path, "getclinical.R"),  local = TRUE)$value
    source(file.path(server.path, "survival.R"),  local = TRUE)$value
    source(file.path(server.path, "volcano.R"),  local = TRUE)$value
    source(file.path(server.path, "heatmap.R"),  local = TRUE)$value
    source(file.path(server.path, "dmr.R"),  local = TRUE)$value
    source(file.path(server.path, "meanMet.R"),  local = TRUE)$value
    source(file.path(server.path, "dea.R"),  local = TRUE)$value
    source(file.path(server.path, "pathview.R"),  local = TRUE)$value
    source(file.path(server.path, "eaplot.R"),  local = TRUE)$value
    source(file.path(server.path, "oncoprint.R"),  local = TRUE)$value
    source(file.path(server.path, "maftools.R"),  local = TRUE)$value
    source(file.path(server.path, "starburst.R"),  local = TRUE)$value
    source(file.path(server.path, "elmer.R"),  local = TRUE)$value
    source(file.path(server.path, "manageSE.R"),  local = TRUE)$value
    source(file.path(server.path, "getinference.R"),  local = TRUE)$value

    # Directory management
    # Config
    # We will create by deafalt a TCGAbiolinksGUI
    dir.create(paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI"), showWarnings = FALSE)
    setwd(file.path(Sys.getenv("HOME"), "TCGAbiolinksGUI"))
    shinyDirChoose(input, 'workingDir',
                   roots=get.volumes(),
                   session=session,
                   restrictions=system.file(package='base'))

    shinyjs::hide("greetbox-outer")

    observe({
        if(!is.null(input$test)) stopApp()  # stop shiny
    })
    # Configuration tab
    output$wd <- renderPrint({
        path <- parseDirPath(get.volumes(isolate({input$workingDir})), input$workingDir)
        if(identical(path, character(0))) path <- getwd()
        return(path)
    })

    observe({
        output$downloadDataBut <- downloadHandler(
            filename <-   function() {
                as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), input$downloadfile)$datapath)
            },
            content <- function(file){
                file.copy(as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), input$downloadfile)$datapath),
                          file)
            }
        )
    })
    #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
    # File selection
    #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
    observe({
        shinyFileChoose(input,
                        'downloadfile',
                        roots = get.volumes(input$workingDir),
                        session = session,
                        restrictions = system.file(package='base'),
                        filetypes = c('csv', 'rda'))
    })

    hide("loading-content", TRUE, "fade")
}
