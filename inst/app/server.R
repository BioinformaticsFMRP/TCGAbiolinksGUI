library(shiny)
library(shinyFiles)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(UpSetR)
library(ggplot2)
library(shinyBS)
library(stringr)
library(ggrepel)
library(pathview)
library(htmlwidgets)
library(ELMER)
library(googleVis)
library(readr)
library(data.table)
library(grid)
options(shiny.maxRequestSize=-1)


getDataCategory <- function(legacy){
    data.category.hamonirzed <- c("Transcriptome Profiling","Copy Number Variation",
                                  "Simple Nucleotide Variation","Simple Nucleotide Variation",
                                  "Raw Sequencing Data","Biospecimen","Clinical")

    data.category.legacy <- c("Copy number variation",
                              "Simple Nucleotide Variation","Simple Nucleotide Variation",
                              "Raw Sequencing Data","Biospecimen","Clinical","Protein expression","Gene expression",
                              "DNA methylation","Raw Microarray Data","Structural Rearrangement","Other")
    if(legacy) return(data.category.legacy)
    return(data.category.hamonirzed)
}

getFileType <-  function(legacy, data.category){
    file.type <- NULL
    if(grepl("Copy number variation",data.category, ignore.case = TRUE) & legacy)
        file.type <- c("nocnv_hg18.seg","nocnv_hg19.seg","hg19.seg","hg18.seg")
    if(grepl("Gene expression", data.category, ignore.case = TRUE)  & legacy)
        file.type <- c("normalized_results","results")

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

getMafTumors <- function(){
    root <- "https://gdc-api.nci.nih.gov/data/"
    maf <- fread("https://gdc-docs.nci.nih.gov/Data/Release_Notes/Manifests/GDC_open_MAFs_manifest.txt",
                 data.table = FALSE, verbose = FALSE, showProgress = FALSE)
    tumor <- unlist(lapply(maf$filename, function(x){unlist(str_split(x,"\\."))[2]}))
    disease <-  gsub("TCGA-","",TCGAbiolinks:::getGDCprojects()$project_id)
    names(disease) <-  TCGAbiolinks:::getGDCprojects()$disease_type
    disease <- sort(disease)
    ret <- disease[disease %in% tumor]
    return(ret)
}

getPlatform <-  function(legacy, data.category){
    platform <- NULL
    if(!legacy) return(platform) # platform is not used for harmonized
    if(grepl("Copy number variation",data.category, ignore.case = TRUE)) platform <- "Affymetrix SNP Array 6.0"
    if(data.category == "Protein expression") platform <- "MDA RPPA Core"
    if(data.category == "Gene expression") platform <- c("Illumina HiSeq","HT_HG-U133A","AgilentG4502A_07_2","AgilentG4502A_07_1","HuEx-1_0-st-v2")
    if(data.category == "DNA methylation") platform <- c("Illumina Human Methylation 450","Illumina Human Methylation 27",
                                                         "Illumina DNA Methylation OMA003 CPI","Illumina DNA Methylation OMA002 CPI",
                                                         "Illumina Hi Seq")
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

#' @title  Server side
#' @description Server side
#' @param input - input signal
#' @param output - output signal
#' @importFrom downloader download
#' @import pathview ELMER TCGAbiolinks SummarizedExperiment shiny ggrepel UpSetR
#' @keywords internal
TCGAbiolinksGUIServer <- function(input, output, session) {
    #addClass(selector = "body", class = "sidebar-collapse")
    setwd(Sys.getenv("HOME"))
    #volumes <- c('Working directory'=getwd())
    #switch(Sys.info()[['sysname']],
    #       Windows= {root = "C:\\"},
    #       Linux  = {root = "~"},
    #       Darwin = {root = file.path("~", "Desktop")})
    dir.create(paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI"), showWarnings = FALSE)

    volumes <- c(TCGAbiolinksGUI=paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI"),home=getwd(),getVolumes()(),temp=tempdir(),wd="./")
    shinyjs::hide("greetbox-outer")

    #-------------------------------------------------------------------------
    #                            TCGA Search
    #-------------------------------------------------------------------------
    #--------------------- START controlling show/hide states ----------------
    observeEvent(input$clinicalIndexed, {
        if(input$clinicalIndexed){
            shinyjs::show("tcgaClinicalTypeFilter")
            shinyjs::hide("tcgaClinicalFilter")
        } else {
            shinyjs::hide("tcgaClinicalTypeFilter")
            shinyjs::show("tcgaClinicalFilter")
        }
    })

    #----------------------- END controlling show/hide states -----------------
    observeEvent(input$tcgaSearchBt, {
        updateCollapse(session, "collapseTCGA", open = "TCGA search results")
        output$tcgaview <- renderGvis({
            closeAlert(session,"tcgaAlert")

            #------------------- STEP 1: Argument-------------------------
            # Set arguments for GDCquery, if value is empty we will set it to FALSE (same as empty)
            tumor <- isolate({input$tcgaProjectFilter})
            data.category <- isolate({input$tcgaDataCategoryFilter})

            # Data type
            data.type <- isolate({input$tcgaDataTypeFilter})
            if(str_length(data.type) == 0) data.type <- FALSE

            # platform
            platform <- isolate({input$tcgaPlatformFilter})
            if(str_length(platform) == 0) platform <- FALSE

            # workflow
            workflow.type <- isolate({input$tcgaWorkFlowFilter})
            if(str_length(workflow.type) == 0) workflow.type <- FALSE

            # file.type
            file.type <- isolate({input$tcgaFileTypeFilter})
            if(str_length(file.type) == 0) file.type <- FALSE

            # access: we will only work with open access data
            #access <- isolate({input$tcgaAcessFilter})
            #if(str_length(access) == 0) access <- FALSE
            access <- "open"

            legacy <- isolate({as.logical(input$tcgaDatabase)})

            # bacode
            text.samples <- isolate({input$tcgaDownloadBarcode})
            if(!is.null(text.samples)){
                barcode <- parse.textarea.input(text.samples)
            }
            if(str_length(barcode) == 0) barcode <- FALSE

            # Samples type
            sample.type <- isolate({input$tcgasamplestypeFilter})

            if(is.null(sample.type)) {
                sample.type <- FALSE
            } else if(str_length(sample.type) == 0) {
                sample.type <- FALSE
            }

            experimental.strategy <- isolate({input$tcgaExpStrategyFilter})
            if(str_length(experimental.strategy) == 0) experimental.strategy <- FALSE

            if(is.null(tumor)){
                createAlert(session, "tcgasearchmessage", "tcgaAlert", title = "Data input error", style =  "danger",
                            content = "Please select a project", append = FALSE)
                return(NULL)
            }
            if(is.null(data.category)){
                createAlert(session, "tcgasearchmessage", "tcgaAlert", title = "Data input error", style =  "danger",
                            content = "Please select a data.category", append = FALSE)
                return(NULL)
            }
            # print(tumor)
            # print(data.category)
            # print(workflow.type)
            # print(legacy)
            # print(platform)
            # print(file.type)
            # print(access)
            # print(barcode)
            # print(experimental.strategy)
            # print(data.type)
            # print(sample.type)

            withProgress(message = 'Search in progress',
                         detail = 'This may take a while...', value = 0, {
                             query = tryCatch({
                                 query <- GDCquery(project = tumor,
                                                   data.category = data.category,
                                                   workflow.type = workflow.type,
                                                   legacy = legacy,
                                                   platform = platform,
                                                   file.type = file.type,
                                                   access = access,
                                                   barcode = barcode,
                                                   experimental.strategy = experimental.strategy,
                                                   sample.type = sample.type,
                                                   data.type = data.type)
                             }, error = function(e) {
                                 createAlert(session, "tcgasearchmessage", "tcgaAlert", title = "Error", style =  "danger",
                                             content = "No results for this query", append = FALSE)
                                 return(NULL)
                             })
                         })
            # print(is.null(query))
            if(is.null(query)) return(NULL)
            not.found <- c()
            tbl <- data.frame()
            results <- query$results[[1]]
            if(any(duplicated(results$cases)))
                createAlert(session, "tcgasearchmessage", "tcgaAlert", title = "Warning", style =  "warning",
                            content = "There are more than one file for the same case.", append = FALSE)

            #------------------- STEP 2: Clinical features-------------------------
            # In this step we will download indexed clinical data
            # Get the barcodes that repects user inputs
            # And select only the results that respects it

            clinical <- GDCquery_clinic(tumor)
            stage <- isolate({input$tcgaClinicalTumorStageFilter})
            stage.idx <- NA
            if(!is.null(stage) & all(str_length(stage) > 0)){
                stage.idx <- sapply(stage, function(y) clinical$tumor_stage %in% y)
                stage.idx <- apply(stage.idx,1,any)
            }

            vital.status <- isolate({input$tcgaClinicalVitalStatusFilter})
            vital.status.idx <- NA
            if(!is.null(vital.status) & all(str_length(vital.status) > 0)){
                vital.status.idx <- sapply(vital.status, function(y) clinical$vital_status %in% y)
                vital.status.idx <- apply(vital.status.idx,1,any)
            }
            race <- isolate({input$tcgaClinicalRaceFilter})
            race.idx <- NA
            if(!is.null(race) & all(str_length(race) > 0)){
                race.idx <- sapply(race, function(y) clinical$race %in% y)
                race.idx <- apply(race.idx,1,any)
            }
            gender <- isolate({input$tcgaClinicalGenderFilter})
            gender.idx <- NA
            if(!is.null(gender) & all(str_length(gender) > 0)){
                gender.idx <- sapply(gender, function(y) clinical$gender %in% y)
                gender.idx <- apply(gender.idx,1,any)
            }
            idx <- apply(data.frame(gender.idx,race.idx,vital.status.idx,stage.idx),1,function(x)all(x,na.rm = TRUE))
            clinical <- clinical[idx,]
            results <- results[substr(results$cases,1,12) %in% clinical$bcr_patient_barcode,]
            # filter samples

            #------------------- STEP 3: Plot-------------------------
            # Plot informations to help the user visualize the data that will
            # be downloaded

            data.type <- gvisPieChart(as.data.frame(table(results$data_type)),
                                      options=list( title="Data type"))
            tissue.definition <- gvisPieChart(as.data.frame(table(results$tissue.definition)),
                                              options=list( title="Tissue definition"))
            experimental_strategy  <- gvisPieChart(as.data.frame(table(results$experimental_strategy )),
                                                   options=list( title="Experimental strategy"))
            analysis.workflow_type <- gvisPieChart(as.data.frame(table(paste0(
                results$analysis$workflow_type,"(", results$analysis$workflow_link, ")"))),
                options=list( title="Workflow type"))
            data.category <- gvisPieChart(as.data.frame(table(results$data_category)),
                                          options=list( title="Data category"))

            gender.plot <- gvisPieChart(as.data.frame(table(clinical$gender)),
                                        options=list( title = "Gender"))
            race.plot <- gvisPieChart(as.data.frame(table(clinical$race)),
                                      options=list( title = "Race"))
            vital.status.plot <- gvisPieChart(as.data.frame(table(clinical$vital_status)),
                                              options=list( title="Vital status"))
            tumor.stage.plot <- gvisPieChart(as.data.frame(table(clinical$tumor_stage)),
                                             options=list( title="Tumor stage"))

            clinical.plots <- gvisMerge(gvisMerge(tumor.stage.plot,vital.status.plot , horizontal = TRUE),
                                        gvisMerge(race.plot, gender.plot, horizontal = TRUE))

            if(!legacy) {
                data.plots <- gvisMerge(gvisMerge(tissue.definition, data.type, horizontal = TRUE),
                                        gvisMerge(experimental_strategy, analysis.workflow_type, horizontal = TRUE))

            } else {
                data.plots <- gvisMerge(gvisMerge(tissue.definition, data.type, horizontal = TRUE),
                                        experimental_strategy)
            }
            gvisMerge(data.plots,clinical.plots)
        })
    })

    observeEvent(input$tcgaPrepareBt,{

        # Dir to save the files
        getPath <- parseDirPath(volumes, input$workingDir)
        if (length(getPath) == 0) getPath <- paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI")

        withProgress(message = 'Download in progress',
                     detail = 'This may take a while...', value = 0, {
                            GDCdownload(query = query, directory = getPath)
                         GDCprepare(query, save = TRUE,
                                    save.filename = ,
                                    summarizedExperiment = isolate({input$prepareRb}),
                                    directory = getPath)
                             createAlert(session, "tcgasearchmessage", "tcgaAlert", title = "Download completed", style =  "success",
                                         content =  paste0("Saved in: ", getPath), append = FALSE)
                             closeAlert(session, "tcgaAlert")
                             createAlert(session, "tcgasearchmessage", "tcgaAlert", title = "Prepare completed", style =  "success",
                                         content =  paste0("Saved in: ", getPath,"<br><ul>", paste(all.saved, collapse = ""),"</ul>"), append = FALSE)

                         })
    })

    # Subtype

    observeEvent(input$tcgaSubtypeBt, {
        updateCollapse(session, "collapseTCGASubtype", open = "Subtype data")
        output$tcgaSubtypetbl <- renderDataTable({
            tumor <- isolate({input$tcgasubtypeFilter})
            tbl <- data.frame()


            result = tryCatch({
                tbl <- rbind(tbl, TCGAquery_subtype(tumor = tumor))


            }, error = function(e) {
                createAlert(session, "tcgaSubtypemessage", "tcgaSubtypeAlert", title = "No results found", style =  "warning",
                            content = "Sorry there are subtypes for your query.", append = FALSE)
            })

            if(is.null(tbl)){
                createAlert(session, "tcgaSubtypemessage", "tcgaSubtypeAlert", title = "No results found", style =  "warning",
                            content = "Sorry there are subtypes for your query.", append = FALSE)
                return()
            } else if(nrow(tbl) ==0) {
                createAlert(session, "tcgaSubtypemessage", "tcgaSubtypeAlert", title = "No results found", style =  "warning",
                            content = "Sorry there are subtypes for your query.", append = FALSE)
                return()
            } else {
                closeAlert(session, "tcgaSubtypeAlert")
                doi <- c("acc"="Comprehensive Pan-Genomic Characterization of Adrenocortical Carcinoma<br>doi:10.1016/j.ccell.2016.04.002",
                         "aml"="Genomic and epigenomic landscapes of adult de novo acute<br>doi:10.1056/NEJMoa1301689",
                         "blca"="Comprehensive molecular characterization of urothelial bladder <br>doi:10.1038/nature12965",
                         "brca"="Comprehensive molecular portraits of human breast tumours<br>doi:10.1038/nature11412",
                         "coad"="Comprehensive molecular characterization of human colon<br>doi:10.1038/nature11252",
                         "gbm"="Molecular Profiling Reveals Biologically Discrete Subsets and Pathways of Progression in Diffuse Glioma.<br>doi:10.1016/j.cell.2015.12.028",
                         "lgg"="Molecular Profiling Reveals Biologically Discrete Subsets and Pathways of Progression in Diffuse Glioma.<br>doi:10.1016/j.cell.2015.12.028",
                         "hnsc"="Comprehensive genomic characterization of head and neck<br>doi:10.1038/nature14129",
                         "kich"="The somatic genomic landscape of chromophobe renal cell carcinoma<br>doi:10.1016/j.ccr.2014.07.014",
                         "kirc"="Comprehensive molecular characterization of clear cell renal cell carcinoma<br>doi:10.1038/nature12222",
                         "kirp"="Comprehensive Molecular Characterization of Papillary Renal-Cell Carcinoma<br>doi:10.1056/NEJMoa1505917",
                         "lihc"="",
                         "luad"="Comprehensive molecular profiling of lung adenocarcinoma <br>doi:10.1038/nature13385",
                         "lusc"="Comprehensive genomic characterization of squamous cell lung cancers<br>doi:10.1038/nature11404",
                         "ovca"= "Integrated genomic analyses of ovarian carcinoma<br>doi:10.1038/nature10166",
                         "pancan"="Multiplatform analysis of 12 cancer types reveals molecular <br>doi:10.1016/j.cell.2014.06.049",
                         "prad"="The Molecular Taxonomy of Primary Prostate Cancer<br>doi:10.1016/j.cell.2015.10.025",
                         "skcm"="Genomic Classification of Cutaneous Melanoma<br>doi:10.1016/j.cell.2015.05.044",
                         "stad"="Comprehensive molecular characterization of gastric adenocarcinoma<br>doi:10.1038/nature13480",
                         "thca"="Integrated genomic characterization of papillary thyroid carcinoma<br>doi:10.1016/j.cell.2014.09.050",
                         "ucec"="Integrated genomic characterization of endometrial carcinoma<br>doi:10.1038/nature12113",
                         "ucs"="")
                if (isolate({input$saveSubtypeRda}) || isolate({input$saveSubtypeCsv})) {
                    save.message <- ""
                    getPath <- parseDirPath(volumes, input$workingDir)
                    if (length(getPath) == 0) getPath <- paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI")
                    filename <- file.path(getPath,paste0(tumor,"_subtype.rda"))
                    if (isolate({input$saveSubtypeRda})) {
                        save(tbl, file = filename)
                        save.message <- paste0(save.message,"<br> File created: ", filename)
                    }
                    if (isolate({input$saveSubtypeCsv})) {
                        write.csv2(tbl, file = gsub("rda","csv",filename))
                        save.message <- paste0(save.message,"<br> File created: ", gsub("rda","csv",filename))
                    }
                    createAlert(session, "tcgaSubtypemessage", "tcgaSubtypeAlert", title = paste0("Success"), style =  "success",
                                content = paste0(save.message,"<br>Source of the data:", doi[tumor]), append = TRUE)

                } else {
                    createAlert(session, "tcgaSubtypemessage", "tcgaSubtypeAlert", title = "Source of the data", style =  "success",
                                content = paste0(doi[tumor]), append = TRUE)
                }

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
        updateCollapse(session, "collapseTCGAClinical", open = "Clinical data table", close = "Information about clinical data")
        output$tcgaClinicaltbl <- renderDataTable({
            closeAlert(session, "tcgaClinicalAlert")
            project <- isolate({input$tcgatumorClinicalFilter})
            type <- isolate({input$tcgaClinicalTypeFilter})
            parser <- isolate({input$tcgaClinicalFilter})
            text.samples <- isolate({input$clinicalBarcode})
            tbl <- data.frame()

            getPath <- parseDirPath(volumes, input$workingDir)
            if (length(getPath) == 0) getPath <- paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI")

            withProgress( message = 'Search in progress',
                          detail = 'This may take a while...', value = 0, {
                              result = tryCatch({
                                  if(isolate({input$clinicalIndexed})){
                                      tbl <- rbind(tbl, GDCquery_clinic(project = project, type = type))
                                  } else {
                                      type <- "Clinical"
                                      query <- GDCquery(project = project, data.category = type)
                                      incProgress(1/2, detail = paste("Downloading data"))
                                      result = tryCatch({
                                          GDCdownload(query, directory = getPath)
                                      } , error = function(e) {
                                          createAlert(session, "tcgaClinicalmessage", "tcgaClinicalAlert", title = "No results found", style =  "error",
                                                      content = "There was a problem to download the data. Please try again later.", append = FALSE)

                                      })
                                      clinical <- GDCprepare_clinic(query,parser, directory = getPath)
                                      tbl <- rbind(tbl, clinical)
                                  }
                              }, error = function(e) {
                                  createAlert(session, "tcgaClinicalmessage", "tcgaClinicalAlert", title = "No results found", style =  "warning",
                                              content = "Sorry, we could not find the clinical information for your query.", append = FALSE)
                              })
                          })
            if(is.null(tbl)){
                createAlert(session, "tcgaClinicalmessage", "tcgaClinicalAlert", title = "No results found", style =  "warning",
                            content = "Sorry there are clinical files for your query.", append = FALSE)
                return()
            } else if(nrow(tbl) == 0) {
                createAlert(session, "tcgaClinicalmessage", "tcgaClinicalAlert", title = "No results found", style =  "warning",
                            content = "Sorry there are clinical files for your query.", append = FALSE)
                return()
            } else {
                closeAlert(session, "tcgaClinicalAlert")
                # Saving data
                if (isolate({input$saveClinicalRda}) || isolate({input$saveClinicalCsv})){
                    if(isolate({input$clinicalIndexed})){
                        filename <- paste0(project,"_",type,".rda")
                    } else {
                        filename <- paste0(project,"_clinical_",parser,".rda")
                    }

                    filename <- file.path(getPath,filename)
                    save.message <- ""
                    if (isolate({input$saveClinicalRda})){
                        save(tbl, file = filename)
                        save.message <-  paste0(save.message,"File created: ",filename,"<br>")
                    }
                    if (isolate({input$saveClinicalCsv})){
                        if(grepl("biospecimen",type, ignore.case = TRUE))  tbl$portions <- NULL
                        write_csv(tbl,gsub("rda","csv",filename))
                        save.message <-  paste0(save.message,"File created: ",gsub("rda","csv",filename))
                    }
                    createAlert(session, "tcgaClinicalmessage", "tcgaClinicalAlert", title = "File created", style =  "info",
                                content = save.message, append = TRUE)
                }
                return(tbl)
            }

        },
        options = list(pageLength = 10,
                       scrollX = TRUE,
                       jQueryUI = TRUE,
                       pagingType = "full",
                       autoWidth = TRUE,
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
    # Maf Search
    observeEvent(input$tcgaMafSearchBt, {

        output$tcgaMutationtbl <- renderDataTable({
            tumor <- isolate({input$tcgaMafTumorFilter})
            getPath <- parseDirPath(volumes, input$workingDir)
            if (length(getPath) == 0) getPath <- paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI")

            withProgress(message = 'Download in progress',
                         detail = 'This may take a while...', value = 0, {
                             tbl <- GDCquery_Maf(tumor, directory = getPath, save.csv =  isolate({input$saveMafcsv}))
                         })
            if(is.null(tbl)){
                createAlert(session, "tcgaMutationmessage", "tcgaMutationAlert", title = "No results found", style =  "warning",
                            content = "Sorry there are no results for your query.", append = FALSE)
                return()
            } else if(nrow(tbl) ==0) {
                createAlert(session, "tcgaMutationmessage", "tcgaMutationAlert", title = "No results found", style =  "warning",
                            content = "Sorry there are no results for your query.", append = FALSE)
                return()
            } else {
                closeAlert(session, "tcgaMutationAlert")
            }
            root <- "https://gdc-api.nci.nih.gov/data/"
            maf <- fread("https://gdc-docs.nci.nih.gov/Data/Release_Notes/Manifests/GDC_open_MAFs_manifest.txt",
                         data.table = FALSE, verbose = FALSE, showProgress = FALSE)
            maf <- maf[grepl(tumor,maf$filename),]
            fout <- file.path(getPath, maf$filename)
            closeAlert(session, "tcgaMutationAlert")
            if( isolate({input$saveMafcsv})) {
                createAlert(session, "tcgaMutationmessage", "tcgaMutationAlert", title = "Download completed", style = "success",
                            content =  paste0("Saved file: <br><ul>",fout,"</ul><ul>", gsub(".gz",".csv",fout), "</ul>"), append = FALSE)
            } else {
                createAlert(session, "tcgaMutationmessage", "tcgaMutationAlert", title = "Download completed", style = "success",
                            content =  paste0("Saved file: ",fout), append = FALSE)
            }
            return(tbl)
        },
        options = list(pageLength = 10,
                       scrollX = TRUE,
                       jQueryUI = TRUE,
                       pagingType = "full",
                       lengthMenu = list(c(6, 20, -1), c('10', '20', 'All')),
                       #autoWidth = TRUE,
                       #columnDefs = list(list(width = '200px', targets = 110:ncol(tbl))),
                       #initComplete = JS(
                       #   "function(settings, json) {",
                       #   "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
                       #   "}"),
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


    #----------------------------------------------------------------------
    #                                         MAF
    #----------------------------------------------------------------------

    #-------------------------START controlling show/hide states -----------------
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
    observeEvent(input$maffile, {
        if(!is.null(mut())){
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

    #-------------------------END controlling show/hide states -----------------

    shinyFileChoose(input, 'maffile', roots=volumes, session=session,
                    restrictions=system.file(package='base'), filetypes=c('', 'maf'))
    shinyFileChoose(input, 'mafAnnotation', roots=volumes, session=session,
                    restrictions=system.file(package='base'), filetypes=c('', 'csv','rda'))
    shinyFileChoose(input, 'oncoGenesFiles', roots=volumes, session=session, restrictions=system.file(package='base'))

    annotation.maf <- function(){
        inFile <- input$mafAnnotation
        if (is.null(inFile)) return(NULL)
        file  <- as.character(parseFilePaths(volumes, input$mafAnnotation)$datapath)
        if(tools::file_ext(file)=="csv"){
            se <- read.csv2(file,header = T,stringsAsFactors = FALSE, row.names = 1)
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
    }
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

    genesByFile <- function(){
        inFile <- input$oncoGenesFiles
        if (is.null(inFile)) return(NULL)
        file  <- as.character(parseFilePaths(volumes, inFile)$datapath)
        if(tools::file_ext(file)=="csv"){
            df <- read.csv2(file,header = T,stringsAsFactors = FALSE)
            rownames(df) <- df[,1]
            df[,1] <- NULL
        } else if(tools::file_ext(file)=="rda"){
            df <- get(load(file))
        } else if(tools::file_ext(file)=="txt"){
            df <- read.table(file,header = T)
        } else {
            createAlert(session, "oncomessage", "oncoAlert", title = "Data input error", style =  "danger",
                        content = paste0("Sorry, but I'm expecting a csv, rda or txt file, but I got a: ",
                                         tools::file_ext(file)), append = FALSE)
            return(NULL)
        }
        genes <- NULL
        # if a data frame return the column with gene symbols
        if(class(df)==class(data.frame())){
            if("mRNA" %in% colnames(df)){
                if("status" %in% colnames(df)) df <- subset(df,df$status != "Insignificant")
                aux <- strsplit(df$mRNA,"\\|")
                genes <- unlist(lapply(aux,function(x) x[1]))
            } else if("Gene_symbol" %in% colnames(df)){
                genes <- df$Gene_symbol
            } else {
                createAlert(session, "oncomessage", "oncoAlert", title = "Data input error", style =  "danger",
                            content = paste0("Sorry, but I'm expecting a column called Gene_symbol "), append = FALSE)
                return(NULL)
            }
        }
        return(genes)
    }

    observeEvent(input$oncoprintPlot , {
        closeAlert(session,"oncoAlert")
        output$oncoploting <- renderPlot({
            mut <- isolate({mut()})
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
            } else if( length(genes) > 200){
                createAlert(session, "oncomessage", "oncoAlert", title = "Errors", style =  "danger",
                            content = paste0("The limit of the genes is 200 \n You gave me: ",length(genes)), append = TRUE)
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
            } else{
                closeAlert(session, "oncoAlert")
            }
            withProgress(message = 'Creating plot',
                         detail = 'This may take a while...', value = 0, {

                             TCGAvisualize_oncoprint(mut=mut,genes=genes,annotation = annotation,
                                                     annotation.position=isolate(input$mafAnnotationpos),
                                                     rm.empty.columns = rm.empty.cols,
                                                     show.column.names = show.col.names,
                                                     show.row.barplot = show.row.barplot,
                                                     label.font.size = isolate(input$oncoTextLabelSize),
                                                     rows.font.size = isolate(input$oncoTextRowSize),
                                                     dist.row =  isolate(input$oncoHSpace),
                                                     dist.col =  isolate(input$oncoWSpace),
                                                     row.order =  isolate(input$oncoRowSort),
                                                     annotation.legend.side = isolate(input$oncoAnnotationLegendSide),
                                                     heatmap.legend.side = isolate(input$oncoHeatmapLegendSide),
                                                     color = c("background" = isolate(input$colBG),
                                                               "SNP"=isolate(input$colSNP),"INS"=isolate(input$colINS),
                                                               "DEL"=isolate(input$colDEL),"DNP"=isolate(input$colDNP)))

                         })
        })})
    observeEvent(input$oncoprintPlot , {
        updateCollapse(session, "collapseOnco", open = "Oncoprint")
        output$oncoPlot <- renderUI({
            plotOutput("oncoploting", width = paste0(isolate({input$oncowidth}), "%"), height = isolate({input$oncoheight}))
        })})
    observe({
        updateSelectizeInput(session, 'mafAnnotationcols', choices = as.character(colnames(annotation.maf())), server = TRUE)
    })
    observe({
        updateSelectizeInput(session, 'tcgaDataCategoryFilter', choices =  getDataCategory(as.logical(input$tcgaDatabase)), server = TRUE)
    })

    observe({
        updateSelectizeInput(session, 'tcgaDataTypeFilter', choices =  getDataType(as.logical(input$tcgaDatabase),input$tcgaDataCategoryFilter), server = TRUE)
        if(is.null(getDataType(as.logical(input$tcgaDatabase),input$tcgaDataCategoryFilter))) {
            shinyjs::hide("tcgaDataTypeFilter")
        } else {
            shinyjs::show("tcgaDataTypeFilter")
        }
    })
    observe({
        updateSelectizeInput(session, 'tcgaPlatformFilter', choices =  getPlatform(as.logical(input$tcgaDatabase),input$tcgaDataCategoryFilter), server = TRUE)
        if(is.null(getPlatform(as.logical(input$tcgaDatabase),input$tcgaDataCategoryFilter))) {
            shinyjs::hide("tcgaPlatformFilter")
        } else {
            shinyjs::show("tcgaPlatformFilter")
        }
    })
    observe({
        updateSelectizeInput(session, 'tcgaWorkFlowFilter', choices =  getWorkFlow(as.logical(input$tcgaDatabase),input$tcgaDataCategoryFilter), server = TRUE)
        if(is.null(getWorkFlow(as.logical(input$tcgaDatabase),input$tcgaDataCategoryFilter))) {
            shinyjs::hide("tcgaWorkFlowFilter")
        } else {
            shinyjs::show("tcgaWorkFlowFilter")
        }
    })

    observe({
        updateSelectizeInput(session, 'tcgaExpStrategyFilter', choices =  getExpStrategy(as.logical(input$tcgaDatabase),input$tcgaPlatformFilter), server = TRUE)
        if(is.null(getExpStrategy(as.logical(input$tcgaDatabase),input$tcgaPlatformFilter))) {
            shinyjs::hide("tcgaExpStrategyFilter")
        } else {
            shinyjs::show("tcgaExpStrategyFilter")
        }
    })

    observe({
        updateSelectizeInput(session, 'tcgaMafTumorFilter', choices =  getMafTumors(), server = TRUE)
    })

    observe({
        updateSelectizeInput(session, 'tcgaFileTypeFilter', choices =  getFileType(as.logical(input$tcgaDatabase),input$tcgaDataCategoryFilter), server = TRUE)
        if(is.null(getFileType(as.logical(input$tcgaDatabase),input$tcgaDataCategoryFilter))) {
            shinyjs::hide("tcgaFileTypeFilter")
        } else {
            shinyjs::show("tcgaFileTypeFilter")
        }
    })


    observe({
        print(input$tcgaProjectFilter)
        tryCatch({
            clin <- GDCquery_clinic(input$tcgaProjectFilter)
            updateSelectizeInput(session, 'tcgaClinicalGenderFilter', choices =  unique(clin$gender), server = TRUE)
            updateSelectizeInput(session, 'tcgaClinicalVitalStatusFilter', choices =  unique(clin$vital_status), server = TRUE)
            updateSelectizeInput(session, 'tcgaClinicalRaceFilter', choices =  unique(clin$race), server = TRUE)
            updateSelectizeInput(session, 'tcgaClinicalTumorStageFilter', choices =  unique(clin$tumor_stage), server = TRUE)
            shinyjs::show("tcgaClinicalGenderFilter")
            shinyjs::show("tcgaClinicalVitalStatusFilter")
            shinyjs::show("tcgaClinicalRaceFilter")
            shinyjs::show("tcgaClinicalTumorStageFilter")
        }, error = function(e){
            shinyjs::hide("tcgaClinicalGenderFilter")
            shinyjs::hide("tcgaClinicalVitalStatusFilter")
            shinyjs::hide("tcgaClinicalRaceFilter")
            shinyjs::hide("tcgaClinicalTumorStageFilter")
        })
    })

    observe({
        updateSelectizeInput(session, 'oncoGenes', choices = as.character(mut()$Hugo_Symbol), server = TRUE)
    })

    ##----------------------------------------------------------------------
    #                             Volcano plot
    ##----------------------------------------------------------------------
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

    shinyFileChoose(input, 'volcanofile', roots=volumes, session=session, restrictions=system.file(package='base'),filetypes=c('excel', 'csv'))

    volcanodata <-  reactive({

        inFile <- input$volcanofile
        if (is.null(inFile)) return(NULL)
        file  <- as.character(parseFilePaths(volumes, inFile)$datapath)
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
                         df <- read.csv2(file,header = T,row.names = 1)
                     })
        return(df)

    })

    observeEvent(input$volcanofile, {
        file  <- basename(as.character(parseFilePaths(volumes, input$volcanofile)$datapath))
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
            file  <- basename(as.character(parseFilePaths(volumes, input$volcanofile)$datapath))
            file <- unlist(str_split(file,"_"))
            groupCol <- file[3]
            group1 <- file[4]
            group2 <- file[5]
            names <- NULL

            # methylation pipeline
            if(isolate({input$volcanoInputRb})=="met"){

                diffcol <- paste("diffmean", group1, group2,sep = ".")
                pcol <- paste("p.value.adj", group1, group2,sep = ".")

                if(isolate({input$volcanoNames})) names <- data$probeID
                label <- c("Not Significant",
                           "Hypermethylated",
                           "Hypomethylated")
                label[2:3] <-  paste(label[2:3], "in", group2)

                # Update data into a file

                statuscol <- paste("status",group1,group2,sep = ".")
                statuscol2 <- paste("status",group2,group1,sep = ".")
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
                    csv <- paste0(paste("DMR_results",
                                        gsub("_",".",groupCol),
                                        gsub("_",".",group1),
                                        gsub("_",".",group2),
                                        "pcut",y.cut,
                                        "meancut",x.cut,
                                        sep = "_"),
                                  ".csv")
                    write.csv2(data,file =  csv)
                    createAlert(session, "volcanomessage", "volcanoAlert", title = "File created", style =  "success",
                                content = paste0(getwd(),"/",csv), append = FALSE)

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
                if(isolate({input$volcanoNames})) names <- as.character(data$mRNA)
                data$Gene_Symbol  <- as.character(data$mRNA)
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
                    out.filename <- paste0(paste("DEA_results",
                                                 gsub("_",".",groupCol),
                                                 gsub("_",".",group1),
                                                 gsub("_",".",group2),
                                                 "pcut", y.cut,
                                                 "logFC.cut",x.cut,
                                                 sep="_"),
                                           ".csv")
                    write.csv2(data, file = out.filename)
                    createAlert(session, "volcanomessage", "volcanoAlert", title = "File created", style =  "success",
                                content =  paste0(getwd(),"/",out.filename), append = FALSE)
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
    # automatically change the type based in the input
    observe({
        if(!is.null(input$volcanofile)){
            file  <- basename(as.character(parseFilePaths(volumes, input$volcanofile)$datapath))
            selected <- "met"
            if(grepl("DEA",file))  selected <- "exp"
            updateRadioButtons(session, "volcanoInputRb", selected = selected)
        }
    })
    observe({
        if(!is.null(input$volcanoHighlight)){
            updateCheckboxInput(session, "volcanoNames",  value = TRUE)
            updateCheckboxInput(session, "volcanoNamesFill",  value = TRUE)
            updateSelectizeInput(session, 'volcanoShowHighlitgh', selected = "highlighted")
        } else {
            updateSelectizeInput(session, 'volcanoShowHighlitgh', selected = "significant")
        }
    })
    observe({
        data <- volcanodata()
        if(!is.null(data)) {
            file  <- basename(as.character(parseFilePaths(volumes, input$volcanofile)$datapath))
            if(grepl("DEA",file)){
                if("Gene_Symbol" %in% colnames(data)){
                    updateSelectizeInput(session, 'volcanoHighlight', choices = as.character(na.omit(unique(data$Gene_Symbol))), server = TRUE)
                }
                if("mRNA" %in% colnames(data)){
                    updateSelectizeInput(session, 'volcanoHighlight', choices = as.character(na.omit(unique(data$mRNA))), server = TRUE)
                }
            } else {
                updateSelectizeInput(session, 'volcanoHighlight', choices = as.character(na.omit(unique(data$probeID))), server = TRUE)
            }
        }
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

    ##----------------------------------------------------------------------
    #                             DMR analysis
    ##----------------------------------------------------------------------

    #-------------------------START controlling show/hide states -----------------
    observeEvent(input$meanmetSetLimits, {
        if(input$meanmetSetLimits){
            shinyjs::show("meanmetylimup")
            shinyjs::show("meanmetylimlow")
        } else {
            shinyjs::hide("meanmetylimup")
            shinyjs::hide("meanmetylimlow")
        }
    })
    observeEvent(input$meanmetSortCB, {
        if(input$meanmetSortCB){
            shinyjs::show("meanmetsort")
        } else {
            shinyjs::hide("meanmetsort")
        }
    })
    #-------------------------END controlling show/hide states -----------------
    observeEvent(input$dmrAnalysis , {

        groups <- t(combn(isolate({input$dmrgroups}),2))
        print(groups)
        # read the data from the downloaded path
        # prepare it
        se <- isolate({dmrdata()})
        se <- subset(se,subset = (rowSums(is.na(assay(se))) == 0))
        withProgress(message = 'DMR analysis in progress',
                     detail = 'This may take a while...', value = 0, {
                         message <- "<br>Saving the results also in a csv file:<ul>"
                         for(i in 1:nrow(groups)) {
                             incProgress(1/(nrow(groups)+ 1 ), detail = paste(groups[i,1]," vs ", groups[i,2]))
                             group1 <- groups[i,1]
                             group2 <- groups[i,2]
                             se <- TCGAanalyze_DMR(data = se,
                                                   groupCol = isolate({input$dmrgroupCol}),
                                                   group1 = group1,
                                                   group2 = group2,
                                                   plot.filename = paste0("DMR_volcano_",group1,"_vs_",group2,".pdf"),
                                                   p.cut = isolate({input$dmrpvalue}),
                                                   diffmean.cut = isolate({input$dmrthrsld}),
                                                   cores = isolate({input$dmrcores}))
                             message <- paste0(message,"<li>DMR_results_",
                                               gsub("_",".",isolate({input$dmrgroupCol})),
                                               "_", gsub("_",".",group1), "_", gsub("_",".",group2), "_",
                                               "pcut_",isolate({input$dmrpvalue}), "_",
                                               "meancut_",isolate({input$dmrthrsld}),".csv</li>")
                         }
                         file  <- as.character(parseFilePaths(volumes, input$dmrfile)$datapath)
                         if(!grepl("results",file)) file <- gsub(".rda","_results.rda",file)
                         save(se,file = file)
                         incProgress(1/(nrow(groups) + 1 ), detail = paste("Saving results"))
                     })
        createAlert(session, "dmrmessage", "dmrAlert", title = "DMR completed", style =  "danger",
                    content = paste0("Summarized Experiment object with results saved in: ", file, message,"<ul>"),
                    append = FALSE)
    })
    shinyFileChoose(input, 'dmrfile', roots=volumes, session=session,
                    restrictions=system.file(package='base'),filetypes=c('', 'rda'))
    shinyFileChoose(input, 'meanmetfile', roots=volumes, session=session,
                    restrictions=system.file(package='base'),filetypes=c('', 'rda'))
    shinyFileChoose(input, 'heatmapfile', roots=volumes, session=session,
                    restrictions=system.file(package='base'),filetypes=c('', 'rda'))
    shinyFileChoose(input, 'heatmapresultsfile', roots=volumes, session=session,
                    restrictions=system.file(package='base'),filetypes=c('', 'csv'))

    meandata <-  reactive({
        inFile <- input$meanmetfile
        if (is.null(inFile)) return(NULL)
        file  <- as.character(parseFilePaths(volumes, input$meanmetfile)$datapath)

        withProgress(message = 'Loading data',
                     detail = 'This may take a while...', value = 0, {
                         result.file <- gsub(".rda","_results.rda",file)
                         if(file.exists(result.file)) {
                             se <- get(load(result.file))
                         } else {
                             se <- get(load(file))
                         }
                     })
        if(class(se)!= class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
            createAlert(session, "dmrmessage", "dmrAlert", title = "Data input error", style =  "danger",
                        content = paste0("Sorry, but I'm expecting a Summarized Experiment object, but I got a: ",
                                         class(se)), append = FALSE)
            return(NULL)
        }
        return(se)

    })
    dmrdata <-  reactive({
        inFile <- input$dmrfile
        if (is.null(inFile)) return(NULL)
        file  <- as.character(parseFilePaths(volumes, input$dmrfile)$datapath)

        withProgress(message = 'Loading data',
                     detail = 'This may take a while...', value = 0, {
                         result.file <- gsub(".rda","_results.rda",file)
                         if(file.exists(result.file)) {
                             se <- get(load(result.file))
                         } else {
                             se <- get(load(file))
                         }
                     })
        if(class(se)!= class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
            createAlert(session, "dmrmessage", "dmrAlert", title = "Data input error", style =  "danger",
                        content = paste0("Sorry, but I'm expecting a Summarized Experiment object, but I got a: ",
                                         class(se)), append = FALSE)
            return(NULL)
        }
        return(se)

    })

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

    observeEvent(input$dmrgroupCol , {
        updateSelectizeInput(session, 'dmrgroups', choices = {
            if (class(dmrdata()) == class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
                if (!is.null(dmrdata()) & input$dmrgroupCol != "" )
                    as.character(colData(dmrdata())[,input$dmrgroupCol])
            }}, server = TRUE)
    })
    observe({
        updateSelectizeInput(session, 'dmrgroupCol', choices = {
            # remove numeric columns
            data <- dmrdata()
            if(!is.null(data)){
                data <- colData(data)
                #data <- data[,which(apply(data,2,function(x) class(parse_guess((x)))) == "character")]
                as.character(colnames(data))
            }
        }, server = TRUE)
    })
    observe({
        updateSelectizeInput(session, 'meanmetsubgroupCol', choices = {
            data <- meandata()
            if(!is.null(data)){
                data <- colData(data)
                #data <- data[,which(apply(data,2,function(x) class(parse_guess((x)))) == "character")]
                as.character(colnames(data))
            }
        }, server = TRUE)
    })
    observe({
        updateSelectizeInput(session, 'meanmetgroupCol', choices = {
            data <- meandata()
            if(!is.null(data)){
                data <- colData(data)
                #data <- data[,which(apply(data,2,function(x) class(parse_guess((x)))) == "character")]
                as.character(colnames(data))
            }
        }, server = TRUE)
    })

    observeEvent(input$meanmetPlot , {

        output$mean.plotting <- renderPlot({
            closeAlert(session, "meanmetAlert")
            jitter <- isolate({input$meanmetplotjitter})
            sort <- NULL
            if(isolate({input$meanmetSortCB})) sort <- isolate({input$meanmetsort})
            angle <- isolate({input$meanmetAxisAngle})
            data <- meandata()

            y.limits <- NULL
            if(isolate({input$meanmetSetLimits})) y.limits <- c(isolate({input$meanmetylimlow}),isolate({input$meanmetylimup}))

            if(is.null(data)){
                createAlert(session, "meanmetmessage", "meanmetAlert", title = "Missing data", style =  "danger",
                            content = paste0("Please select the data"), append = FALSE)
                return(NULL)
            }

            if(isolate({input$meanmetgroupCol}) == "") {
                group <- NULL
            } else {
                group <- isolate({input$meanmetgroupCol})
            }

            if(is.null(group)){
                createAlert(session, "meanmetmessage", "meanmetAlert", title = "Missing group column", style =  "danger",
                            content = paste0("Please select group column"), append = FALSE)
                return(NULL)
            }

            if(isolate({input$meanmetsubgroupCol}) == "") {
                subgroup <- NULL
            } else {
                subgroup <- isolate({input$meanmetsubgroupCol})
            }

            legend.position <- isolate({input$meanmetLegendPos})
            legend.title.position <- isolate({input$meanmetLegendTitlePos})
            legend.ncols <- isolate({input$meanmetLegendncols})
            add.axis.x.text <- isolate({input$meanmetXnames})
            withProgress(message = 'Creating plot',
                         detail = 'This may take a while...', value = 0, {
                             if(is.null(sort)){
                                 TCGAvisualize_meanMethylation(data=data,
                                                               groupCol=group,
                                                               subgroupCol=subgroup,
                                                               filename = NULL,
                                                               y.limits = y.limits,
                                                               legend.position = legend.position,
                                                               legend.title.position = legend.title.position,
                                                               legend.ncols = legend.ncols,
                                                               add.axis.x.text = add.axis.x.text,
                                                               plot.jitter = jitter,
                                                               axis.text.x.angle = angle )
                             } else {
                                 TCGAvisualize_meanMethylation(data=data,
                                                               groupCol=group,
                                                               subgroupCol=subgroup,
                                                               filename = NULL,
                                                               legend.position = legend.position,
                                                               legend.title.position = legend.title.position,
                                                               legend.ncols = legend.ncols,
                                                               add.axis.x.text = add.axis.x.text,
                                                               y.limits = y.limits,
                                                               plot.jitter = jitter,
                                                               axis.text.x.angle = angle,
                                                               sort=sort)
                             }
                         })
        })})

    observeEvent(input$meanmetPlot , {
        updateCollapse(session, "collapsemeanmet", open = "Mean DNA methylation plot")
        output$meanMetplot <- renderUI({
            plotOutput("mean.plotting", width = paste0(isolate({input$meanmetwidth}), "%"), height = isolate({input$meanmetheight}))
        })})

    output$probesSE <- renderDataTable({
        data <- dmrdata()

        if(!is.null(data)) {
            df <- as.data.frame(values(data))

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

    ##----------------------------------------------------------------------
    #                             Heatmap
    ##----------------------------------------------------------------------

    #-------------------------START controlling show/hide states -----------------

    observe({
        if(input$heatmapTypeInputRb == "met") {
            shinyjs::show("heatmapProbesInputRb")
            shinyjs::hide("heatmapGenesInputRb")
            shinyjs::hide("heatmapGenesTextArea")
            shinyjs::hide("heatmap.upGenesCb")
            shinyjs::hide("heatmap.downGenewsCb")
            if(input$heatmapProbesInputRb == "text"){
                shinyjs::show("heatmapProbesTextArea")
                shinyjs::hide("heatmap.hypoprobesCb")
                shinyjs::hide("heatmap.hyperprobesCb")
            } else {
                shinyjs::hide("heatmapProbesTextArea")
                shinyjs::show("heatmap.hypoprobesCb")
                shinyjs::show("heatmap.hyperprobesCb")
            }
        } else if(input$heatmapTypeInputRb == "exp") {
            shinyjs::show("heatmapGenesInputRb")
            shinyjs::hide("heatmapProbesTextArea")
            shinyjs::hide("heatmap.hypoprobesCb")
            shinyjs::hide("heatmap.hyperprobesCb")
            shinyjs::hide("heatmapProbesInputRb")
            if(input$heatmapGenesInputRb == "text"){
                shinyjs::show("heatmapGenesTextArea")
                shinyjs::hide("heatmap.upGenesCb")
                shinyjs::hide("heatmap.downGenewsCb")
            } else {
                shinyjs::hide("heatmapGenesTextArea")
                shinyjs::show("heatmap.upGenesCb")
                shinyjs::show("heatmap.downGenewsCb")
            }
        }
    })
    observeEvent(input$heatmap.sortCb, {
        shinyjs::toggle("heatmapSortCol")
    })

    #-------------------------END controlling show/hide states -----------------

    heatmapresultdata <-  reactive({
        inFile <- input$heatmapresultsfile
        if (is.null(inFile)) return(NULL)
        file  <- as.character(parseFilePaths(volumes, inFile)$datapath)
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
                         df <- read.csv2(file,header = TRUE, stringsAsFactors = FALSE)
                     })
        return(df)
    })

    heatmapdata <-  reactive({
        inFile <- input$heatmapfile
        if (is.null(inFile)) return(NULL)
        file  <- as.character(parseFilePaths(volumes, input$heatmapfile)$datapath)

        withProgress(message = 'Loading data',
                     detail = 'This may take a while...', value = 0, {
                         result.file <- gsub(".rda","_results.rda",file)
                         if(file.exists(result.file)) {
                             se <- get(load(result.file))
                         } else {
                             se <- get(load(file))
                         }
                     })
        if(class(se)!= class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
            createAlert(session, "heatmapmessage", "heatmapAlert", title = "Data input error", style =  "danger",
                        content = paste0("Sorry, but I'm expecting a Summarized Experiment object, but I got a: ",
                                         class(se)), append = FALSE)
            return(NULL)
        }
        return(se)

    })
    observeEvent(input$heatmapPlotBt , {
        output$heatmap.plotting <- renderPlot({

            # get information from file
            file  <- basename(as.character(parseFilePaths(volumes, input$heatmapresultsfile)$datapath))
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

            if( isolate({input$heatmapTypeInputRb}) == "met") type <-  "methylation"
            if( isolate({input$heatmapTypeInputRb}) == "exp") type <-  "expression"

            if(nchar(sortCol) ==0 &  isolate({input$heatmap.sortCb})){
                createAlert(session, "heatmapmessage", "heatmapAlert", title = "Columns metadata", style =  "danger",
                            content = paste0("Please select the heatmapSortCol"),append = FALSE)
                return(NULL)
            }

            if( isolate({input$heatmapTypeInputRb})=="met"){
                # ---------------- probes selection
                if(isolate({input$heatmapProbesInputRb}) == "Status"){
                    sig.probes <- ""
                    if(isolate({input$heatmap.hypoprobesCb})) sig.probes <- c("Hypomethylated")
                    if(isolate({input$heatmap.hyperprobesCb})) sig.probes <- c("Hypermethylated",sig.probes)
                    sig.probes <- paste(sig.probes,"in",group2)
                    # Get hypo methylated and hypermethylated probes
                    idx <- paste("status",group1,group2, sep=".")
                    print(table(results.data[,idx]))
                    probes <- results.data[,idx] %in% sig.probes
                } else {
                    sig.probes <- parse.textarea.input(isolate({input$heatmapProbesTextArea}))
                    probes <- which(results.data$probeID %in% sig.probes)
                }
                data <- data[probes,]
                results.data <- results.data[probes,]
            } else {
                if(isolate({input$heatmapGenesInputRb}) == "Status"){
                    sig.genes <- ""
                    if(isolate({input$heatmap.upGenesCb})) sig.genes <- c("Upregulated")
                    if(isolate({input$heatmap.downGenewsCb})) sig.genes <- c("Downregulated",sig.genes)
                    sig.genes <- paste(sig.genes,"in",group2)
                    # Get hypo methylated and hypermethylated probes
                    genes <- results.data[,"status"] %in% sig.genes
                } else {
                    sig.genes <- parse.textarea.input(isolate({input$heatmapGenesTextArea}))
                    aux <- strsplit(results.data$mRNA,"\\|")
                    results.data$gene <- unlist(lapply(aux,function(x) x[2]))
                    genes <- which(results.data$gene %in% sig.genes)
                }

                data <- data[genes,]
                results.data <- results.data[genes,]

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
                if(length(rowmdata) > 0) row.metadata <- subset(results.data, select=c(rowmdata))
            }
            withProgress(message = 'Creating plot',
                         detail = 'This may take a while...', value = 0, {
                             if(!isolate({input$heatmap.sortCb})) {
                                 p <-  TCGAvisualize_Heatmap(data=assay(data),
                                                             col.metadata=col.metadata,
                                                             row.metadata=row.metadata,
                                                             title = "Heatmap",
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
                                                             title = "Heatmap",
                                                             cluster_rows = cluster_rows,
                                                             show_column_names = show_column_names,
                                                             cluster_columns = cluster_columns,
                                                             show_row_names = show_row_names,
                                                             sortCol = sortCol,
                                                             type = type,
                                                             scale = scale)
                             }
                             incProgress(1/2)
                             ComplexHeatmap::draw(p)
                         })
        })})

    observeEvent(input$heatmapPlotBt , {
        updateCollapse(session, "collapseHeatmap", open = "Heatmap")
        output$heatmapPlot <- renderUI({
            plotOutput("heatmap.plotting", width = paste0(isolate({input$heatmapwidth}), "%"), height = isolate({input$heatmapheight}))
        })})

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
    })

    #--------------------------------------------------------------------
    #                                    EA
    #--------------------------------------------------------------------
    #--------------------- START controlling show/hide states ----------------
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
    #----------------------- END controlling show/hide states -----------------
    shinyFileChoose(input, 'eaGenesFiles', roots=volumes, session=session, restrictions=system.file(package='base'))
    eaGenesByFile <- function(){
        inFile <- input$eaGenesFiles
        if (is.null(inFile)) return(NULL)
        file  <- as.character(parseFilePaths(volumes, inFile)$datapath)
        if(tools::file_ext(file)=="csv"){
            df <- read.csv2(file,header = T,stringsAsFactors = FALSE)
            rownames(df) <- df[,1]
            df[,1] <- NULL
        } else if(tools::file_ext(file)=="rda"){
            df <- get(load(file))
        } else if(tools::file_ext(file)=="txt"){
            df <- read.table(file,header = T)
        } else {
            createAlert(session, "oncomessage", "oncoAlert", title = "Data input error", style =  "danger",
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
                createAlert(session, "oncomessage", "oncoAlert", title = "Data input error", style =  "danger",
                            content = paste0("Sorry, but I'm expecting a column called Gene_symbol "), append = FALSE)
                return(NULL)
            }
        }
        return(genes)
    }
    observeEvent(input$eaplot , {
        updateCollapse(session, "collapseEA", open = "EA plots")
        output$eaPlot <- renderUI({
            plotOutput("ea.plotting", width = paste0(isolate({input$eawidth}), "%"), height = isolate({input$eaheight}))
        })})

    observeEvent(input$eaplot , {
        output$ea.plotting <- renderPlot({

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
                    createAlert(session, "oncomessage", "oncoAlert", title = "Data input error", style =  "danger",
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

                             nbPlot <- table(c(is.null(ResMF),is.null(ResBP),is.null(ResCC),is.null(ResPat)))["FALSE"]

                             if(nbPlot == 1)  mfrow = c(1,1)
                             if(nbPlot == 2)  mfrow = c(2,1)
                             if(nbPlot > 2)  mfrow = c(2,2)
                             if(is.na(nbPlot)) return (NULL)
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
                                                     nRGTab = genes,
                                                     text.size = isolate({input$eaSizeText}),
                                                     xlim = xlim,
                                                     mfrow = mfrow,
                                                     nBar = isolate({input$nBar}),
                                                     filename = NULL)
                         })
        })})
    #--------------------------------------------------------------------------
    #                           Profile plot
    #--------------------------------------------------------------------------
    #--------------------- START controlling show/hide states -----------------
    shinyjs::hide("profileplotgroup")
    shinyjs::hide("profileplotsubtype")
    shinyjs::hide("profileplotrmnagroup")
    shinyjs::hide("profileplotrmnasub")
    observeEvent(input$profileplotfile, {
        if(!is.null(profileplotdata())){
            shinyjs::show("profileplotgroup")
            shinyjs::show("profileplotsubtype")
            shinyjs::show("profileplotrmnagroup")
            shinyjs::show("profileplotrmnasub")
        }
    })
    #----------------------- END controlling show/hide states -----------------

    observe({
        data <- profileplotdata()
        if(!is.null(data)) {
            # remove numeric columns
            #data <- data[,which(apply(data,2,function(x) class(parse_guess((x)))) == "character")]
            names <- colnames(data)[apply(data,2,function(x) length(unique(x))) > 1]
            updateSelectizeInput(session, 'profileplotgroup', choices = {
                as.character(names)
            }, server = TRUE)
            updateSelectizeInput(session, 'profileplotsubtype', choices = {
                as.character(names)
            }, server = TRUE)
        }
    })


    profileplotdata <- function(){
        inFile <- input$profileplotfile
        if (is.null(inFile)) return(NULL)
        file  <- as.character(parseFilePaths(volumes, inFile)$datapath)
        if(tools::file_ext(file)=="csv"){
            df <- read.csv2(file,header = T, row.names = 1)
        } else if(tools::file_ext(file)=="rda"){
            df <- get(load(file))
        } else {
            createAlert(session, "profileplotmessage", "profileplotAlert", title = "Data input error", style =  "danger",
                        content = paste0("Sorry, but I'm expecting a csv or rda file, but I got a: ",
                                         tools::file_ext(file)), append = FALSE)
            return(NULL)
        }

        if(class(df)!= class(data.frame())){
            createAlert(session, "profileplotmessage", "profileplotAlert", title = "Data input error", style =  "danger",
                        content = paste0("Sorry, but I'm expecting a Data frame object, but I got a: ",
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
        if (x == 4) { m1 <- -5.0; m3 <- -1.5}
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
            closeAlert(session, "profileplotAlert")

            data <- isolate({profileplotdata()})
            subtypeCol <- isolate({input$profileplotsubtype})
            groupCol <-  isolate({input$profileplotgroup})
            na.rm.groups <-  isolate({input$profileplotrmnagroup})
            na.rm.subtypes  <-  isolate({input$profileplotrmnasub})
            m1  <-  isolate({input$margin1})
            m2  <-  isolate({input$margin2})
            m3  <-  isolate({input$margin3})
            m4  <-  isolate({input$margin4})

            axis.title.size <-  isolate({input$profileplot.axis.title.size})
            axis.textsize <-  isolate({input$profileplot.axis.textsize})
            legend.size <-  isolate({input$profileplot.legend.size})
            legend.title.size <- isolate({input$profileplot.legend.title.size})
            geom.label.size <- isolate({input$profileplot.geom.label.size})

            if(is.null(data)){
                createAlert(session, "profileplotmessage", "profileplotAlert", title = "Missing data", style =  "danger",
                            content = paste0("Please select the data"), append = FALSE)
                return(NULL)
            }

            if(is.null(groupCol) || nchar(groupCol) == 0){
                createAlert(session, "profileplotmessage", "profileplotAlert", title = "Missing group selection", style =  "danger",
                            content = paste0("Please select the group column"), append = FALSE)
                return(NULL)
            }

            if(is.null(subtypeCol) || nchar(subtypeCol) == 0){
                createAlert(session, "profileplotmessage", "profileplotAlert", title = "Missing subgroup selection", style =  "danger",
                            content = paste0("Please select the subgroup column"), append = FALSE)
                return(NULL)
            }


            withProgress(message = 'Creating plot',
                         detail = 'This may take a while...', value = 0, {

                             TCGAvisualize_profilePlot(data = data,
                                                       groupCol=groupCol,
                                                       subtypeCol=subtypeCol,
                                                       na.rm.groups = na.rm.groups,
                                                       na.rm.subtypes = na.rm.subtypes,
                                                       plot.margin=c(m1,m2,m3,m4),
                                                       axis.title.size=axis.title.size,
                                                       axis.textsize=axis.textsize,
                                                       legend.size=legend.size,
                                                       legend.title.size=legend.title.size,
                                                       geom.label.size = geom.label.size,
                                                       geom.label.color = isolate({input$profileplotColorTextLeftBar}))

                         })
        })})

    observeEvent(input$profileplotBt , {
        updateCollapse(session, "collapseprofileplot", open = "Profile plot")
        output$profileplot <- renderUI({
            plotOutput("profile.plotting", width = paste0(isolate({input$profilewidth}), "%"), height = isolate({input$profileheight}))
        })})
    shinyFileChoose(input, 'profileplotfile', roots=volumes, session=session,
                    restrictions=system.file(package='base'),filetypes=c('', 'rda','csv'))


    #------------------------------------------------
    # Survival plot
    # -----------------------------------------------
    #--------------------- START controlling show/hide states -----------------
    #shinyjs::hide("survivalplotgroup")
    #shinyjs::hide("survivalplotMain")
    #shinyjs::hide("survivalplotLegend")
    #shinyjs::hide("survivalplotLimit")
    #shinyjs::hide("survivalplotPvalue")
    #observeEvent(input$survivalplotfile, {
    #    if(!is.null(survivalplotdata())){
    #        shinyjs::show("survivalplotgroup")
    #        shinyjs::show("survivalplotMain")
    #        shinyjs::show("survivalplotLegend")
    #        shinyjs::show("survivalplotLimit")
    #        shinyjs::show("survivalplotPvalue")
    #    }
    #})
    #----------------------- END controlling show/hide states -----------------
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
        closeAlert(session, "survivalAlert")
        inFile <- input$survivalplotfile
        if (is.null(inFile)) return(NULL)
        file  <- as.character(parseFilePaths(volumes, inFile)$datapath)
        if(tools::file_ext(file)=="csv"){
            df <- read.csv2(file,header = T)
            rownames(df) <- df[,1]
            df[,1] <- NULL
        } else if(tools::file_ext(file)=="rda"){
            df <- get(load(file))
        } else {
            createAlert(session, "survivalmessage", "survivalAlert", title = "Data input error", style =  "danger",
                        content = paste0("Sorry, but I'm expecting a csv or rda file, but I got a: ",
                                         tools::file_ext(file)), append = FALSE)
            return(NULL)
        }
        if(class(df)!= class(data.frame())){
            createAlert(session, "survivalmessage", "survivalAlert", title = "Data input error", style =  "danger",
                        content = paste0("Sorry, but I'm expecting a Data frame object, but I got a: ",
                                         class(df)), append = FALSE)
            return(NULL)
        }
        cols <- c("days_to_death","days_to_last_followup","vital_status")
        if(!(all(cols %in% colnames(df)))){
            createAlert(session, "survivalmessage", "survivalAlert", title = "Data input error", style =  "danger",
                        content = paste0("Sorry, but I'm expecting columns: ",paste(cols,collapse = ", ")), append = FALSE)
            return(NULL)
        }
        return(df)
    }
    observeEvent(input$survivalplotgroup , {
        if(isolate({input$survivalplotgroup}) != ""){
            updateTextInput(session, "survivalplotLegend", value = isolate({input$survivalplotgroup}))
        }
    })

    observeEvent(input$survivalplotBt , {
        output$survival.plotting <- renderPlot({

            closeAlert(session, "survivalAlert")

            data <- isolate({survivalplotdata()})
            legend <- isolate({input$survivalplotLegend})
            main <- isolate({input$survivalplotMain})
            clusterCol <-  isolate({input$survivalplotgroup})
            cut.off <- isolate({input$survivalplotLimit})
            print.pvalue <- isolate({input$survivalplotPvalue})


            #---------------------------
            # Input verification
            #---------------------------
            if(is.null(data)){
                createAlert(session, "survivalmessage", "survivalAlert", title = "Missing data", style =  "danger",
                            content = paste0("Please select the data"), append = FALSE)
                return(NULL)
            }

            if(is.null(clusterCol) || nchar(clusterCol) == 0){
                createAlert(session, "survivalmessage", "survivalAlert", title = "Missing group", style =  "danger",
                            content = paste0("Please select group column"), append = FALSE)
                return(NULL)
            }

            if(length(unique(data[,clusterCol])) == 1){
                createAlert(session, "survivalmessage", "survivalAlert", title = "Data input error", style =  "danger",
                            content = paste0("Sorry, but I'm expecting at least two groups"), append = FALSE)
                return(NULL)
            }
            #-=-=-=-=-=-=-=-=--==-=-=-=-=-=-=-=-=

            withProgress(message = 'Creating plot',
                         detail = 'This may take a while...', value = 0, {

                             TCGAanalyze_survival(data = data,
                                                  clusterCol = clusterCol,
                                                  filename = NULL,
                                                  legend = legend,
                                                  main = main,
                                                  cutoff = cut.off,
                                                  print.value = print.pvalue)

                         })
        })})

    observeEvent(input$survivalplotBt , {
        updateCollapse(session, "collapsesurvivalplot", open = "survival plot")
        output$survivalplot <- renderUI({
            plotOutput("survival.plotting", width = paste0(isolate({input$survivalwidth}), "%"), height = isolate({input$survivalheight}))
        })})
    shinyFileChoose(input, 'survivalplotfile', roots=volumes, session=session,
                    restrictions=system.file(package='base'),filetypes=c('', 'csv','rda'))

    # -------------------------------------------------
    # DEA
    # -------------------------------------------------
    #--------------------- START controlling show/hide states -----------------
    #shinyjs::hide("deanormalizationmet")
    #shinyjs::hide("deanormalizationmet")
    observeEvent(input$deanormalization, {
        shinyjs::toggle("deanormalizationmet")
    })
    observeEvent(input$deafilter, {
        shinyjs::toggle("deafilteringmet")
        shinyjs::toggle("deafilteringcut")
    })
    #----------------------- END controlling show/hide states -----------------
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
        withProgress(message = 'Differential Expression Analysis in progress',
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

                         exp <- TCGAanalyze_DEA(mat1 = dataFilt[,samples.g1],
                                                mat2 = dataFilt[,samples.g2],
                                                Cond1type = g1 ,
                                                Cond2type = g2,
                                                #fdr.cut  = fdr.cut,
                                                #logFC.cut = logFC.cut,
                                                method = method)
                         exp <- TCGAanalyze_LevelTab(exp,
                                                     typeCond1 = g1,
                                                     typeCond2 = g2,
                                                     TableCond1 = dataFilt[,samples.g1],
                                                     TableCond2 = dataFilt[,samples.g2])
                         exp$status <- "Insignificant"
                         exp[exp$logFC >= logFC.cut & exp$FDR <= fdr.cut,"status"] <- paste0("Upregulated in ", g2)
                         exp[exp$logFC <= -logFC.cut & exp$FDR <= fdr.cut,"status"] <- paste0("Downregulated in ", g2)
                         exp$Gene_Symbol <- unlist(lapply(strsplit(exp$mRNA,"\\|"),function(x) x[2]))
                     })

        out.filename <- paste0(paste("DEA_results",gsub("_",".",groupCol),
                                     gsub("_",".",g1), gsub("_",".",g2),
                                     "pcut",fdr.cut,"logFC.cut",logFC.cut,sep="_"),".csv")
        write.csv2(exp, file = out.filename)
        createAlert(session, "deamessage", "deaAlert", title = "DEA completed", style =  "danger",
                    content = out.filename, append = FALSE)
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
    pathway.data <- function(){
        inFile <- input$pathewayexpfile
        if (is.null(inFile)) return(NULL)
        file  <- as.character(parseFilePaths(volumes, input$pathewayexpfile)$datapath)
        if(tools::file_ext(file)=="csv"){
            se <- read.csv2(file,header = TRUE, stringsAsFactors = FALSE)
        } else if(tools::file_ext(file)=="rda"){
            se <- get(load(file))
        }

        return(se)
    }
    shinyFileChoose(input, 'pathewayexpfile', roots=volumes, session=session, restrictions=system.file(package='base'))

    observeEvent(input$pathwaygraphBt , {

        data <- pathway.data()
        pathway.id <- isolate({input$pathway.id})
        kegg.native <- isolate({input$kegg.native.checkbt})

        print(head(data))
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

        createAlert(session, "deamessage", "deaAlert", title = "Pathway graph created", style =  "success",
                    content = paste0("Results saved in: ", pathway.id,extension), append = FALSE)

    })

    #------------------------------------------------
    # Starburst
    # -----------------------------------------------
    observeEvent(input$starburstNames, {
        toggle("starburstNamesFill")
    })

    # get DMR result name and update the mean cut and pcut
    observeEvent(input$starburstmetfile, {
        file  <- basename(as.character(parseFilePaths(volumes, input$starburstmetfile)$datapath))
        if(length(file) > 0){
            file <- unlist(str_split(file,"_"))
            group1 <- file[4]
            group2 <- file[5]
            pcut <- file[7]
            meancut <- gsub(".csv","",file[9])
            updateNumericInput(session, "starburstmetdiff", value = meancut)
            updateNumericInput(session, "starburstmetFDR", value = pcut)
        }
    })

    # get DEA result name and update the mean cut and pcut
    observeEvent(input$starburstexpfile, {
        file  <- basename(as.character(parseFilePaths(volumes, input$starburstexpfile)$datapath))
        if(length(file) > 0){
            file <- unlist(str_split(file,"_"))
            group1 <- file[4]
            group2 <- file[5]
            pcut <- file[7]
            fccut <- gsub(".csv","",file[9])
            updateNumericInput(session, "starburstexpFC", value = fccut)
            updateNumericInput(session, "starburstexFDR", value = pcut)
        }
    })

    # Main function
    starburst <- function(){

        closeAlert(session,"starburstAlert")
        if(is.null(result.dea.data())){
            createAlert(session, "starburstmessage", "starburstAlert", title = "Missing data", style =  "danger",
                        content = paste0("Please select the differential expression results"), append = FALSE)
            return(NULL)
        }
        if(is.null(result.dmr.data())){
            createAlert(session, "starburstmessage", "starburstAlert", title = "Missing data", style =  "danger",
                        content = paste0("Please select the differential DNA methylation results"), append = FALSE)
            return(NULL)
        }
        logFC.cut <- isolate({input$starburstexpFC})
        exp.p.cut <- isolate({input$starburstexFDR})
        diffmean.cut <- isolate({input$starburstmetdiff})
        met.p.cut <- isolate({input$starburstmetFDR})
        exp <- result.dea.data()
        met <-result.dmr.data()
        names <- isolate({input$starburstNames})
        names.fill <- isolate({input$starburstNamesFill})
        colors <- c(isolate({input$sbcolInsignigicant}),
                    isolate({input$sbcolUpHypo}),
                    isolate({input$sbcolDownHypo}),
                    isolate({input$sbcolHypo}),
                    isolate({input$sbcolHyper}),
                    isolate({input$sbcolUp}),
                    isolate({input$sbcolDown}),
                    isolate({input$sbcolUpHyper}),
                    isolate({input$sbcolDownHyper}))


        file  <- basename(as.character(parseFilePaths(volumes, isolate({input$starburstmetfile}))$datapath))
        if(length(file) > 0){
            file <- unlist(str_split(file,"_"))
            group1 <- file[4]
            group2 <- file[5]
        }
        file  <- basename(as.character(parseFilePaths(volumes, isolate({input$starburstexpfile}))$datapath))
        if(length(file) > 0){
            file <- unlist(str_split(file,"_"))
            exp.group1 <- file[4]
            exp.group2 <- file[5]
        }

        if(group1 == exp.group1 & group2 == exp.group2){

            result <- TCGAvisualize_starburst(met = met,
                                              exp = exp,
                                              group1 = group1,
                                              group2 = group2,
                                              color = colors,
                                              names = names,
                                              names.fill = names.fill,
                                              exp.p.cut = exp.p.cut,
                                              met.p.cut = met.p.cut,
                                              diffmean.cut = diffmean.cut,
                                              circle = isolate({input$starburstCircle}),
                                              logFC.cut = logFC.cut,
                                              return.plot = TRUE)
            out.filename <- paste0(paste("Starburst_results", group1, group2,
                                         "exp.p.cut", exp.p.cut, "logFC.cut", logFC.cut,
                                         "met.diffmean", diffmean.cut, "met.p.cut", met.p.cut,
                                         sep = "_"),".csv")
            if(isolate({input$starburstSave})) {
                write.csv2(result$starburst, file = out.filename)
                createAlert(session, "starburstmessage", "starburstAlert", title = "Results saved", style =  "success",
                            content = paste0("Results saved in: ", out.filename), append = FALSE)
            }
            return(result)
        }
    }
    # -------------- Starburst plot
    observeEvent(input$starburstPlot , {
        # validate input
        output$starburst.plot <- renderPlot({
            withProgress(message = 'Creating plot',
                         detail = 'This may take a while...', value = 0, {
                             aux <- starburst()
                             if(!is.null(aux)) {
                                 return(aux$plot)
                             }
                         })
        })})

    observeEvent(input$starburstPlot , {
        updateCollapse(session, "collapsedea", open = "dea plots")
        output$starburstPlot <- renderUI({
            plotOutput("starburst.plot", width = paste0(isolate({input$starburstwidth}), "%"), height = isolate({input$starburstheight}))
        })})

    # Starburst plot input data
    result.dea.data <- function(){
        inFile <- input$starburstexpfile
        if (is.null(inFile)) return(NULL)
        file  <- as.character(parseFilePaths(volumes, input$starburstexpfile)$datapath)
        if(tools::file_ext(file)=="csv"){
            se <- read.csv2(file,header = T,row.names = 1)
        } else if(tools::file_ext(file)=="rda"){
            se <- get(load(file))
        }
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
        if(tools::file_ext(file)=="csv"){
            se <- read.csv2(file,header = T, row.names = 1)
        } else if(tools::file_ext(file)=="rda"){
            se <- get(load(file))
        }

        #if(class(se)!= class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
        #    createAlert(session, "deamessage", "deaAlert", title = "Data input error", style =  "danger",
        #                content = paste0("Sorry, but I'm expecting a Summarized Experiment object, but I got a: ",
        #                                 class(se)), append = FALSE)
        #    return(NULL)
        #}
        return(se)
    }
    shinyFileChoose(input, 'starburstmetfile', roots=volumes, session=session,
                    restrictions=system.file(package='base'), filetypes=c('excel', 'csv'))
    shinyFileChoose(input, 'starburstexpfile', roots=volumes, session=session,
                    restrictions=system.file(package='base'), filetypes=c('excel', 'csv'))

    output$starburstResult <- renderDataTable({

        data <- starburst()
        if(!is.null(data)) return(as.data.frame(data$starburst))
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

    #------------------------------------------------
    # ELMER
    # -----------------------------------------------
    #--------------------- START controlling show/hide states -----------------
    #shinyjs::hide("survivalplotgroup")
    #shinyjs::hide("survivalplotMain")
    #shinyjs::hide("survivalplotLegend")
    #shinyjs::hide("survivalplotLimit")
    #shinyjs::hide("survivalplotPvalue")
    observeEvent(input$scatter.plot.type, {
        type <- isolate({input$elmerPlotType})
        if(type =="scatter.plot"){
            scatter.type <- isolate({input$scatter.plot.type})
            if(scatter.type == "tf"){
                shinyjs::show("scatter.plot.tf")
                shinyjs::show("scatter.plot.motif")
                shinyjs::hide("scatter.plot.genes")
                shinyjs::hide("scatter.plot.probes")
                shinyjs::hide("scatter.plot.nb.genes")
            } else if(scatter.type == "pair"){
                shinyjs::hide("scatter.plot.tf")
                shinyjs::hide("scatter.plot.motif")
                shinyjs::show("scatter.plot.genes")
                shinyjs::show("scatter.plot.probes")
                shinyjs::hide("scatter.plot.nb.genes")
            } else {
                shinyjs::hide("scatter.plot.tf")
                shinyjs::hide("scatter.plot.motif")
                shinyjs::hide("scatter.plot.genes")
                shinyjs::show("scatter.plot.probes")
                shinyjs::show("scatter.plot.nb.genes")
            }
        }
    })
    observeEvent(input$elmerPlotType, {
        type <- isolate({input$elmerPlotType})
        if(type =="scatter.plot"){
            shinyjs::show("scatter.plot.type")
            shinyjs::hide("schematic.plot.type")
            shinyjs::hide("schematic.plot.genes")
            shinyjs::hide("schematic.plot.probes")
            shinyjs::hide("ranking.plot.motif")
            shinyjs::hide("ranking.plot.tf")
            scatter.type <- isolate({input$scatter.plot.type})
            if(scatter.type == "tf"){
                shinyjs::show("scatter.plot.tf")
                shinyjs::show("scatter.plot.motif")
                shinyjs::hide("scatter.plot.genes")
                shinyjs::hide("scatter.plot.probes")
                shinyjs::hide("scatter.plot.nb.genes")
            } else if(scatter.type == "pair"){
                shinyjs::hide("scatter.plot.tf")
                shinyjs::hide("scatter.plot.motif")
                shinyjs::show("scatter.plot.genes")
                shinyjs::show("scatter.plot.probes")
                shinyjs::hide("scatter.plot.nb.genes")
            } else {
                shinyjs::hide("scatter.plot.tf")
                shinyjs::hide("scatter.plot.motif")
                shinyjs::hide("scatter.plot.genes")
                shinyjs::show("scatter.plot.probes")
                shinyjs::show("scatter.plot.nb.genes")
            }
        } else if(type =="schematic.plot"){
            shinyjs::hide("scatter.plot.type")
            shinyjs::hide("scatter.plot.tf")
            shinyjs::hide("scatter.plot.motif")
            shinyjs::hide("scatter.plot.genes")
            shinyjs::hide("scatter.plot.probes")
            shinyjs::hide("scatter.plot.nb.genes")
            shinyjs::show("schematic.plot.type")
            shinyjs::hide("ranking.plot.motif")
            shinyjs::hide("ranking.plot.tf")
            type <- isolate({input$schematic.plot.type})
            if(type =="genes"){
                shinyjs::hide("schematic.plot.probes")
                shinyjs::show("schematic.plot.genes")
            } else {
                shinyjs::show("schematic.plot.probes")
                shinyjs::hide("schematic.plot.genes")
            }
        } else if(type =="ranking.plot"){
            shinyjs::hide("scatter.plot.type")
            shinyjs::hide("scatter.plot.tf")
            shinyjs::hide("scatter.plot.motif")
            shinyjs::hide("scatter.plot.genes")
            shinyjs::hide("scatter.plot.probes")
            shinyjs::hide("scatter.plot.nb.genes")
            shinyjs::hide("schematic.plot.type")
            shinyjs::hide("schematic.plot.genes")
            shinyjs::hide("schematic.plot.probes")
            shinyjs::show("ranking.plot.motif")
            shinyjs::show("ranking.plot.tf")
        } else if(type =="motif.enrichment.plot"){
            shinyjs::hide("scatter.plot.type")
            shinyjs::hide("scatter.plot.tf")
            shinyjs::hide("scatter.plot.motif")
            shinyjs::hide("scatter.plot.genes")
            shinyjs::hide("scatter.plot.probes")
            shinyjs::hide("scatter.plot.nb.genes")
            shinyjs::hide("schematic.plot.type")
            shinyjs::hide("schematic.plot.genes")
            shinyjs::hide("schematic.plot.probes")
            shinyjs::hide("ranking.plot.motif")
            shinyjs::hide("ranking.plot.tf")
        }
    })

    observeEvent(input$schematic.plot.type, {
        type <- isolate({input$elmerPlotType})
        if(type =="schematic.plot"){
            type <- isolate({input$schematic.plot.type})
            if(type =="genes"){
                shinyjs::hide("schematic.plot.probes")
                shinyjs::show("schematic.plot.genes")
            } else {
                shinyjs::show("schematic.plot.probes")
                shinyjs::hide("schematic.plot.genes")
            }
        }
    })
    #----------------------- END controlling show/hide states -----------------
    shinyDirChoose(input, 'elmerFolder', roots=volumes, session=session, restrictions=system.file(package='base'))
    shinyFileChoose(input, 'elmermeefile', roots=volumes, session=session,
                    restrictions=system.file(package='base'), filetypes=c('', 'rda'))
    shinyFileChoose(input, 'elmerresultsfile', roots=volumes, session=session,
                    restrictions=system.file(package='base'), filetypes=c('', 'rda'))
    shinyFileChoose(input, 'elmermetfile', roots=volumes, session=session,
                    restrictions=system.file(package='base'), filetypes=c('', 'rda'))
    shinyFileChoose(input, 'elmerexpfile', roots=volumes, session=session,
                    restrictions=system.file(package='base'), filetypes=c('', 'rda'))

    # mee create
    observe({
        data.met <- mee.met()
        updateSelectizeInput(session, 'elmermeetype', choices = {
            if(!is.null(data.met)) as.character(colnames(colData(data.met)))
        }, server = TRUE)
    })
    observeEvent(input$elmermeetype, {
        data.met <- mee.met()
        updateSelectizeInput(session, 'elmermeesubtype', choices = {
            if(!is.null(data.met)) as.character(unique(colData(data.met)[,input$elmermeetype]))
        }, server = TRUE)
    })
    observeEvent(input$elmermeetype, {
        data.met <- mee.met()
        updateSelectizeInput(session, 'elmermeesubtype2', choices = {
            if(!is.null(data.met)) as.character(unique(colData(data.met)[,input$elmermeetype]))
        }, server = TRUE)
    })

    mee.exp <-  reactive({
        inFile <- input$elmerexpfile
        if (is.null(inFile)) return(NULL)
        file  <- as.character(parseFilePaths(volumes, inFile)$datapath)

        withProgress(message = 'Loading data',
                     detail = 'This may take a while...', value = 0, {
                         exp <- get(load(file))
                     })
        return(exp)
    })
    mee.met <-  reactive({
        inFile <- input$elmermetfile
        if (is.null(inFile)) return(NULL)
        file  <- as.character(parseFilePaths(volumes, inFile)$datapath)

        withProgress(message = 'Loading data',
                     detail = 'This may take a while...', value = 0, {
                         met <- get(load(file))
                     })

        if(class(met)!= class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
            createAlert(session, "elmermessage", "elmerAlert", title = "Data input error", style =  "danger",
                        content = paste0("Sorry, but I'm expecting a Summarized Experiment object, but I got a: ",
                                         class(met)), append = FALSE)
            return(NULL)
        }

        return(met)
    })


    observeEvent(input$elmerpreparemee, {
        exp.elmer <- mee.exp()
        met.elmer <- mee.met()
        column <- isolate({input$elmermeetype})
        subtype1 <-  isolate({input$elmermeesubtype})
        subtype2 <-  isolate({input$elmermeesubtype2})

        if(is.null(column)){
            closeAlert(session, "elmerAlert")
            createAlert(session, "elmermessage", "elmerAlert", title = "Type missing", style =  "danger",
                        content =   "Please select the column with type", append = TRUE)
            return(NULL)
        }
        if(is.null(subtype1) | is.null(subtype2)){
            closeAlert(session, "elmerAlert")
            createAlert(session, "elmermessage", "elmerAlert", title = "Subtype missing", style =  "danger",
                        content =   "Please select the two subtypes", append = TRUE)
            return(NULL)
        }
        if(is.null(exp.elmer) | is.null(met.elmer)){
            closeAlert(session, "elmerAlert")
            createAlert(session, "elmermessage", "elmerAlert", title = "Subtype missing", style =  "danger",
                        content =   "Please upload the two summarized Experiment objects", append = TRUE)
            return(NULL)
        }


        # Data: get only samples of subtype1 and subtype2
        exp.elmer <- subset(exp.elmer, select=(colData(exp.elmer)[,column] %in% c(subtype1,subtype2)))
        met.elmer <- subset(met.elmer, select=(colData(met.elmer)[,column] %in%  c(subtype1,subtype2)))

        # Get barcodes for subtype1
        sample.info <-  colData(exp.elmer)
        samples <- sample.info[sample.info[,column] == isolate({input$elmermeesubtype2}),]$barcode
        withProgress(message = 'Creating mee data',
                     detail = 'This may take a while...', value = 0, {

                         exp.elmer <- TCGAprepare_elmer(exp.elmer, platform = "IlluminaHiSeq_RNASeqV2",save = FALSE)
                         incProgress(1/5, detail = paste0('Expression preparation is done'))

                         met.elmer <- TCGAprepare_elmer(met.elmer, platform = "HumanMethylation450",met.na.cut = isolate({input$elmermetnacut}), save = FALSE)
                         incProgress(1/5, detail = paste0('Methylation preparation is done'))

                         geneAnnot <- txs()
                         geneAnnot$GENEID <- paste0("ID",geneAnnot$GENEID)
                         geneInfo <- promoters(geneAnnot,upstream = 0, downstream = 0)
                         probe <- get.feature.probe()

                         # create mee object, use @ to access the matrices inside the object
                         mee <- fetch.mee(meth = met.elmer, exp = exp.elmer, TCGA = TRUE, probeInfo = probe, geneInfo = geneInfo)
                         incProgress(1/5, detail = paste0('Mee is done'))

                         # Relabel samples in the mee object: subtype1 is control
                         mee@sample$TN[mee@sample$ID %in% substr(samples,1,15)] <- "Control"
                         save(mee,file = paste0("mee_",column,"_",gsub(" ","-",subtype1),"_",gsub(" ","-",subtype2),".rda"))
                         incProgress(2/5, detail = paste0('Saving is done'))
                         createAlert(session, "elmermessage", "elmerAlert", title = "Mee created", style =  "success",
                                     content =   paste0("Mee file created: mee_",column,"_",gsub(" ","-",subtype1),"_",gsub(" ","-",subtype2),".rda"), append = TRUE)
                     })
    })
    # Input data
    elmer.results.data <-  reactive({
        inFile <- input$elmerresultsfile
        if (is.null(inFile)) return(NULL)
        file  <- as.character(parseFilePaths(volumes, inFile)$datapath)

        withProgress(message = 'Loading data',
                     detail = 'This may take a while...', value = 0, {
                         load(file,envir = globalenv())
                     })
    })

    meedata <-  reactive({
        inFile <- input$elmermeefile
        if (is.null(inFile)) return(NULL)
        file  <- as.character(parseFilePaths(volumes, inFile)$datapath)

        withProgress(message = 'Loading data',
                     detail = 'This may take a while...', value = 0, {
                         mee <- get(load(file))
                     })
        return(mee)
    })

    # Updates based on uploaded data
    observe({
        updateSelectizeInput(session, 'scatter.plot.probes', choices = {
            if(!is.null(elmer.results.data())) as.character(Sig.probes$probe)
        }, server = TRUE)
    })
    observe({
        updateSelectizeInput(session, 'scatter.plot.tf', choices = {
            if(!is.null(elmer.results.data())) as.character(rownames(TF.meth.cor))
        }, server = TRUE)
    })
    observe({
        updateSelectizeInput(session, 'scatter.plot.motif', choices = {
            if(!is.null(elmer.results.data())) as.character(names(enriched.motif))
        }, server = TRUE)
    })

    observeEvent(input$scatter.plot.probes, {
        updateSelectizeInput(session, 'scatter.plot.genes', choices = {
            if(!is.null(elmer.results.data())) as.character(nearGenes[[input$scatter.plot.probes]]$GeneID)
        }, server = TRUE)
    })

    observe({
        updateSelectizeInput(session, 'ranking.plot.tf', choices = {
            if(!is.null(elmer.results.data())) as.character(rownames(TF.meth.cor))
        }, server = TRUE)
    })

    observe({
        updateSelectizeInput(session, 'ranking.plot.motif', choices = {
            if(!is.null(elmer.results.data())) as.character(colnames(TF.meth.cor))
        }, server = TRUE)
    })

    observe({
        updateSelectizeInput(session, 'schematic.plot.probes', choices = {
            mee <- meedata()
            if(!is.null(elmer.results.data()) & !is.null(mee)){
                pair.obj <- fetch.pair(pair=pair,
                                       probeInfo = getProbeInfo(mee),
                                       geneInfo = getGeneInfo(mee))
                as.character(pair.obj@pairInfo$Probe)
            }
        }, server = TRUE)
    })
    observe({
        updateSelectizeInput(session, 'schematic.plot.genes', choices = {
            mee <- meedata()
            if(!is.null(elmer.results.data()) & !is.null(mee)){
                pair.obj <- fetch.pair(pair=pair,
                                       probeInfo = getProbeInfo(mee),
                                       geneInfo = getGeneInfo(mee))

                as.character(pair.obj@pairInfo$GeneID)
            }
        }, server = TRUE)
    })
    observeEvent(input$elmerAnalysisBt, {
        getPath <- parseDirPath(volumes, isolate({input$elmerFolder}))
        mee <- meedata()
        if(is.null(mee)){
            closeAlert(session, "elmerAlert")
            createAlert(session, "elmermessage", "elmerAlert", title = "Mee object missing", style =  "success",
                        content =   "Please upload the mee object for this plot", append = TRUE)
            return(NULL)
        }
        library(parallel)
        direction <- c("hyper","hypo")

        for (j in direction){
            withProgress(message = 'ELMER analysis',
                         detail = paste0('Direction: ',j), value = 0, {
                             print(j)
                             dir.out <- paste0(getPath,"/elmer/",j)
                             dir.create(dir.out, recursive = TRUE)
                             #--------------------------------------
                             # STEP 3: Analysis                     |
                             #--------------------------------------
                             # Step 3.1: Get diff methylated probes |
                             #--------------------------------------
                             setProgress(value = which(j == direction) * 0.0, message = "Step 1", detail = paste0(j,": get.diff.meth"), session = getDefaultReactiveDomain())
                             Sig.probes <- get.diff.meth(mee,
                                                         cores=isolate({input$elmercores}),
                                                         dir.out =dir.out,
                                                         diff.dir=j,
                                                         pvalue = isolate({input$elmermetpvalue}))

                             #-------------------------------------------------------------
                             # Step 3.2: Identify significant probe-gene pairs            |
                             #-------------------------------------------------------------
                             # Collect nearby 20 genes for Sig.probes
                             setProgress(value = which(j == direction) * 0.1, message = "Step 2", detail = paste0(j,": GetNearGenes"), session = getDefaultReactiveDomain())
                             nearGenes <- GetNearGenes(TRange=getProbeInfo(mee, probe=Sig.probes$probe),
                                                       cores=isolate({input$elmercores}),
                                                       geneAnnot=getGeneInfo(mee),
                                                       geneNum = isolate({input$elmergetpairNumGenes}))


                             pair <- get.pair(mee=mee,
                                              probes=Sig.probes$probe,
                                              nearGenes=nearGenes,
                                              permu.dir=paste0(dir.out,"/permu"),
                                              dir.out=dir.out,
                                              cores=isolate({input$elmercores}),
                                              label= j,
                                              permu.size=isolate({input$elmergetpairpermu}),
                                              Pe = isolate({input$elmergetpairpvalue}),
                                              percentage =  isolate({input$elmergetpairpercentage}),
                                              portion = isolate({input$elmergetpairportion}),
                                              diffExp = isolate({input$elmergetpairdiffExp}))

                             setProgress(value = which(j == direction) * 0.2, message = "Step 3", detail = paste0(j,": fetch.pair"), session = getDefaultReactiveDomain())
                             Sig.probes.paired <- fetch.pair(pair=pair,
                                                             probeInfo = getProbeInfo(mee),
                                                             geneInfo = getGeneInfo(mee))

                             Sig.probes.paired <- read.csv(paste0(dir.out,"/getPair.",j,".pairs.significant.csv"),
                                                           stringsAsFactors=FALSE)[,1]

                             #-------------------------------------------------------------
                             # Step 3.3: Motif enrichment analysis on the selected probes |
                             #-------------------------------------------------------------
                             if(length(Sig.probes.paired) > 0 ){
                                 #-------------------------------------------------------------
                                 # Step 3.3: Motif enrichment analysis on the selected probes |
                                 #-------------------------------------------------------------
                                 setProgress(value = which(j == direction) * 0.3, message = "Step 4", detail = paste0(j,": get.enriched.motif"), session = getDefaultReactiveDomain())
                                 probe <- get.feature.probe()
                                 enriched.motif <- get.enriched.motif(probes = Sig.probes.paired,
                                                                      dir.out = dir.out,
                                                                      label = j,
                                                                      background.probes = probe$name,
                                                                      lower.OR =  isolate({input$elmergetenrichedmotifLoweOR}),
                                                                      min.incidence = isolate({input$elmergetenrichedmotifMinIncidence}))
                                 motif.enrichment <- read.csv(paste0(dir.out,"/getMotif.",j,".motif.enrichment.csv"),
                                                              stringsAsFactors=FALSE)
                                 if(length(enriched.motif) > 0){
                                     #-------------------------------------------------------------
                                     # Step 3.4: Identifying regulatory TFs                        |
                                     #-------------------------------------------------------------
                                     print("get.TFs")
                                     setProgress(value = which(j == direction) * 0.4, message = "Step 5", detail = paste0(j,": get.TFs"), session = getDefaultReactiveDomain())

                                     TF <- get.TFs(mee = mee,
                                                   enriched.motif = enriched.motif,
                                                   dir.out = dir.out,
                                                   cores = isolate({input$elmercores}),
                                                   label = j,
                                                   percentage = isolate({input$elmergetTFpercentage}))
                                     TF.meth.cor <- get(load(paste0(dir.out,"/getTF.",j,".TFs.with.motif.pvalue.rda")))
                                     save(TF, enriched.motif, Sig.probes.paired,
                                          pair, nearGenes, Sig.probes, motif.enrichment, TF.meth.cor,
                                          file=paste0(dir.out,"/ELMER_results_",j,".rda"))
                                 }
                             }
                             setProgress(value = which(j == direction) * 0.5, message = "Step 5", detail = paste0('Analysis in direction completed: ',j), session = getDefaultReactiveDomain())
                         })
        }
    })

    observeEvent(input$elmerPlotBt , {
        output$elmer.plot <- renderPlot({
            plot.type <- isolate({input$elmerPlotType})
            mee <- meedata()
            closeAlert(session, "elmerAlert")

            # Three types:
            # 1 - TF expression vs average DNA methylation
            # 2 - Generate a scatter plot for one probe-gene pair
            # 3 - Generate scatter plots for one probes’ nearby 20 gene expression
            #     vs DNA methylation at this probe
            if(plot.type == "scatter.plot"){
                if(is.null(mee)){
                    createAlert(session, "elmermessage", "elmerAlert", title = "Mee object missing", style =  "success",
                                content =   "Please upload the mee object for this plot", append = TRUE)
                    return(NULL)
                }
                # case 1
                plot.by <- isolate({input$scatter.plot.type})
                if(plot.by == "tf"){
                    if(is.null(isolate({input$scatter.plot.tf}))){
                        closeAlert(session, "elmerAlert")
                        createAlert(session, "elmermessage", "elmerAlert", title = "TFs missing", style =  "success",
                                    content = "Please select two TF", append = TRUE)
                        return(NULL)
                    }

                    if(nchar(isolate({input$scatter.plot.tf})) == 0 | length(isolate({input$scatter.plot.tf})) < 2){
                        closeAlert(session, "elmerAlert")
                        createAlert(session, "elmermessage", "elmerAlert", title = "TFs missing", style =  "success",
                                    content =   "Please select two TF", append = TRUE)
                        return(NULL)
                    }
                    if(nchar(isolate({input$scatter.plot.motif})) == 0){
                        closeAlert(session, "elmerAlert")
                        createAlert(session, "elmermessage", "elmerAlert", title = "Motif missing", style =  "success",
                                    content =   "Please select a motif", append = TRUE)
                        return(NULL)
                    }
                    scatter.plot(mee,byTF=list(TF=isolate({input$scatter.plot.tf}),
                                               probe=enriched.motif[[isolate({input$scatter.plot.motif})]]), category="TN",
                                 save=FALSE,lm_line=TRUE)
                } else if(plot.by == "pair") {
                    if(nchar(isolate({input$scatter.plot.probes})) == 0){
                        closeAlert(session, "elmerAlert")
                        createAlert(session, "elmermessage", "elmerAlert", title = "Probe missing", style =  "success",
                                    content =   "Please select a probe", append = TRUE)
                        return(NULL)
                    }
                    if(nchar(isolate({input$scatter.plot.genes})) == 0){
                        closeAlert(session, "elmerAlert")
                        createAlert(session, "elmermessage", "elmerAlert", title = "Gene missing", style =  "success",
                                    content =   "Please select a gene", append = TRUE)
                        return(NULL)
                    }

                    # case 2
                    scatter.plot(mee,byPair=list(probe=isolate({input$scatter.plot.probes}),gene=c(isolate({input$scatter.plot.genes}))),
                                 category="TN", save=FALSE,lm_line=TRUE)
                } else {
                    # case 3
                    if(nchar(isolate({input$scatter.plot.probes})) == 0){
                        createAlert(session, "elmermessage", "elmerAlert", title = "Probe missing", style =  "success",
                                    content =   "Please select a probe", append = TRUE)
                        return(NULL)
                    }
                    scatter.plot(mee,byProbe=list(probe=isolate({input$scatter.plot.probes}),geneNum=isolate({input$scatter.plot.nb.genes})),
                                 category="TN", dir.out ="./ELMER.example/Result/LUSC", save=FALSE)
                }
            } else if (plot.type == "schematic.plot") {
                if(is.null(mee)){
                    createAlert(session, "elmermessage", "elmerAlert", title = "Mee object missing", style =  "success",
                                content =   "Please upload the mee object for this plot", append = TRUE)
                    return(NULL)
                }
                # Two cases
                # 1 - By probe
                pair.obj <- fetch.pair(pair=pair,
                                       probeInfo = getProbeInfo(mee),
                                       geneInfo = getGeneInfo(mee))
                if(isolate({input$schematic.plot.type}) == "probes"){
                    if(nchar(isolate({input$schematic.plot.probes})) == 0){
                        createAlert(session, "elmermessage", "elmerAlert", title = "Probe missing", style =  "success",
                                    content =   "Please select a probe", append = TRUE)
                        return(NULL)
                    }

                    schematic.plot(pair=pair.obj, byProbe=isolate({input$schematic.plot.probes}),save=FALSE)
                } else if(isolate({input$schematic.plot.type}) == "genes"){
                    if(nchar(isolate({input$schematic.plot.genes})) == 0){
                        createAlert(session, "elmermessage", "elmerAlert", title = "Gene missing", style =  "success",
                                    content =   "Please select a gene", append = TRUE)
                        return(NULL)
                    }

                    # 2 - By genes
                    schematic.plot(pair=pair.obj, byGene=isolate({input$schematic.plot.genes}),save=FALSE)
                }
            } else if(plot.type == "motif.enrichment.plot") {
                motif.enrichment.plot(motif.enrichment=motif.enrichment,
                                      #significant=list(OR=1.3,lowerOR=1.3),
                                      save=FALSE)
            } else if(plot.type == "ranking.plot"){
                if(nchar(isolate({input$ranking.plot.motif})) == 0){
                    createAlert(session, "elmermessage", "elmerAlert", title = "Motif missing", style =  "success",
                                content =   "Please select a motif", append = TRUE)
                    return(NULL)
                }
                label <- list(isolate({input$ranking.plot.tf}))
                names(label) <- isolate({input$ranking.plot.motif})
                gg <- TF.rank.plot(motif.pvalue=TF.meth.cor,
                                   motif=isolate({input$ranking.plot.motif}),
                                   TF.label=label,
                                   save=FALSE)
                # names were not fitting in the plot. Reducing the size
                pushViewport(viewport(height=0.8,width=0.8))
                grid.draw(gg[[1]])
            }
        })
    })
    observeEvent(input$elmerPlotBt , {
        updateCollapse(session, "collapelmer", open = "Plots")
        output$elmerPlot <- renderUI({
            plotOutput("elmer.plot", width = paste0(isolate({input$elmerwidth}), "%"), height = isolate({input$elmerheight}))
        })})

    # Table
    observeEvent(input$elmerTableType , {
        updateCollapse(session, "collapelmer", open = "Results table")
        output$elmerResult <- renderDataTable({
            if(!is.null(elmer.results.data())){
                if(input$elmerTableType == "tf"){
                    as.data.frame(TF)
                } else if(input$elmerTableType == "sigprobes"){
                    as.data.frame(Sig.probes)
                } else if(input$elmerTableType == "motif"){
                    as.data.frame(motif.enrichment)
                } else if(input$elmerTableType == "pair"){
                    as.data.frame(pair)
                }
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
    })

    # Config
    shinyDirChoose(input, 'workingDir', roots=volumes, session=session, restrictions=system.file(package='base'))

    output$wd <- renderPrint({
        path <- parseDirPath(volumes, input$workingDir)
        if (identical(path, character(0))) path <- paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI")
        return(path)
    })
}
