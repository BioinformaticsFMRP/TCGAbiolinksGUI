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
library(ELMER)
library(grid)
options(shiny.maxRequestSize=10000*1024^2)

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
#' @description Server side - Download data from roadmap project
#' @param input - input signal
#' @param output - output signal
#' @importFrom downloader download
#' @keywords internal
biOMICsServer <- function(input, output, session) {
    addClass(selector = "body", class = "sidebar-collapse")
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

        save.dir <- parseDirPath(volumes, input$tcgafolder)
        if(length(save.dir) == 0) {
            filename <- isolate({input$tcgafilename})
        } else {
            filename <- file.path(save.dir,isolate({input$tcgafilename}))
        }

        withProgress( message = 'Prepare in progress',
                      detail = 'This may take a while...', value = 0, {
                          if(!is.null(input$allRows)) {
                              df <- data.frame(name = input$allRows[seq(6, length(input$allRows), 7)])
                              x <- TCGAquery()
                              x <- x[x$name %in% df$name,]

                          } else {

                              tumor <- isolate({input$tcgaTumorFilter})
                              platform <- isolate({input$tcgaExpFilter})
                              level <- isolate({input$tcgaLevelFilter})

                              x <- data.frame()
                              if(length(level) > 0){
                                  for(i in level){
                                      x <- rbind(x,
                                                 TCGAquery(tumor = tumor,
                                                           platform = platform,
                                                           level = i))
                                  }
                              } else {
                                  x <- TCGAquery(tumor = tumor,
                                                 platform = platform,
                                                 level = level)
                              }
                          }
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

    #-------------------------------------------------------------------------
    #                            TCGA Search
    #-------------------------------------------------------------------------
    #--------------------- START controlling show/hide states ----------------
    shinyjs::hide("tcgatumorClinicalFilter")
    shinyjs::hide("tcgaFrnaseqtypeFilter")
    shinyjs::hide("tcgaFrnaseqv2typeFilter")
    shinyjs::hide("tcgaFgwstypeFilter")
    shinyjs::show("addSubTypeTCGA")

    observeEvent(input$tcgaDownloadTypeRb, {
        if(input$tcgaDownloadTypeRb == "none") {
            shinyjs::hide("tcgasamplestypeFilter")
            shinyjs::hide("tcgaDownloadBarcode")
        } else if(input$tcgaDownloadTypeRb == "barcode") {
            shinyjs::hide("tcgasamplestypeFilter")
            shinyjs::show("tcgaDownloadBarcode")
        } else{
            shinyjs::show("tcgasamplestypeFilter")
            shinyjs::hide("tcgaDownloadBarcode")
        }
    })

    observeEvent(input$prepareRb, {
        if(input$prepareRb) {
            shinyjs::show("addSubTypeTCGA")
        } else {
            shinyjs::hide("addSubTypeTCGA")
        }
    })
    observeEvent(input$clinicalSearchType, {
        shinyjs::toggle("clinicalBarcode")
        shinyjs::toggle("tcgatumorClinicalFilter")
    })
    observeEvent(input$tcgaExpFilter, {
        exp <- isolate({input$tcgaExpFilter})

        if(grepl("RNASeqV2", exp,ignore.case = TRUE)) {
            shinyjs::show("tcgaFrnaseqv2typeFilter")
        } else {
            shinyjs::hide("tcgaFrnaseqv2typeFilter")
        }

        if("IlluminaHiSeq_RNASeq" %in% exp) {
            shinyjs::show("tcgaFrnaseqtypeFilter")
        } else {
            shinyjs::hide("tcgaFrnaseqtypeFilter")
        }
        if("Genome_Wide_SNP_6" %in% exp) {
            shinyjs::show("tcgaFgwstypeFilter")
        } else {
            shinyjs::hide("tcgaFgwstypeFilter")
        }
    })
    #observeEvent(input$deafilter, {
    #    shinyjs::toggle("tcgatumorClinicalFilter")
    #    shinyjs::toggle("clinicalBarcode")
    #})
    #----------------------- END controlling show/hide states -----------------
    # Table render
    observeEvent(input$tcgaSearchBt, {
        updateCollapse(session, "collapseTCGA", open = "TCGA search results")
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

    observeEvent(input$tcgafolder , {
        closeAlert(session, "tcgadownloaddirAlert")
        createAlert(session, "tcgaddirmessage", "tcgadownloaddirAlert", title = "Download folder", style =  "success",
                    content =  parseDirPath(volumes, input$tcgafolder), append = TRUE)
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
                         # Decision, if no rows selected download all
                         if(!is.null(input$allRows)) {
                             df <- data.frame(name = input$allRows[seq(6, length(input$allRows), 7)])
                             x <- TCGAquery()
                             x <- x[x$name %in% df$name,]

                         } else {

                             tumor <- isolate({input$tcgaTumorFilter})
                             platform <- isolate({input$tcgaExpFilter})
                             level <- isolate({input$tcgaLevelFilter})

                             x <- data.frame()
                             if(length(level) > 0){
                                 for(i in level){
                                     x <- rbind(x,
                                                TCGAquery(tumor = tumor,
                                                          platform = platform,
                                                          level = i))
                                 }
                             } else {
                                 x <- TCGAquery(tumor = tumor,
                                                platform = platform,
                                                level = level)
                             }
                         }

                         # handling samples input
                         # there are two inputs: type, list of barcodes
                         # if types is selected ignore barcodes,
                         # if barcodes is not null ignore selected
                         text.samples <- isolate({input$tcgaDownloadBarcode})
                         samples <- NULL
                         downloadType <- isolate({input$tcgaDownloadTypeRb})
                         if(downloadType == "barcode"){

                             if(!is.null(text.samples)){
                                 samples <- parse.textarea.input(text.samples)
                             }

                         } else if(downloadType == "type"){
                             samples <- unlist(lapply(samplesType,function(type){
                                 s <- unlist(str_split(aux$barcode,","))
                                 s[grep(type,substr(s,14,15))]
                             }))
                         }

                         if(!is.null(samples)){
                             # filter query
                             idx <- unlist(lapply(samples,function(y) {grep(y,x$barcode)}))
                             x <- x[idx,]
                         }


                         for (i in 1:nrow(x)) {
                             incProgress(1/nrow(x))
                             aux <- x[i,]
                             ftype <- NULL
                             if (grepl("RNASeqV2",aux$Platform,ignore.case = T)) ftype <- rnaseqv2Ftype
                             if (aux$Platform == "IlluminaHiSeq_RNASeq") ftype <- rnaseqFtype
                             if (aux$Platform == "Genome_Wide_SNP_6") ftype <- gwsFtype

                             TCGAdownload(aux, path = getPath,type = ftype,samples = samples)
                         }})
        createAlert(session, "tcgasearchmessage", "tcgaAlert", title = "Download completed", style =  "success",
                    content =  paste0("Saved in: ", getPath), append = FALSE)
    })

    # Subtype

    observeEvent(input$tcgaSubtypeBt, {
        updateCollapse(session, "collapseTCGA", open = "TCGA search results")
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
                if (isolate({input$saveSubtypeRda}) || isolate({input$saveSubtypeCsv})) {
                    save.message <- ""
                    getPath <- parseDirPath(volumes, input$tcgafolder)
                    if (length(getPath) == 0) getPath <- "."
                    filename <- file.path(getPath,paste0(tumor,"_subtype.rda"))
                    if (isolate({input$saveSubtypeRda})) {
                        save(tbl, file = filename)
                        save.message <- paste0(save.message,"<br> File created: ", filename)
                    }
                    if (isolate({input$saveSubtypeCsv})) {
                        write.csv2(tbl, file = gsub("rda","csv",filename))
                        save.message <- paste0(save.message,"<br> File created: ", gsub("rda","csv",filename))
                    }
                    createAlert(session, "tcgasearchmessage", "tcgasearchAlert", title = paste0("Success"), style =  "success",
                                content = paste0(save.message,"<br>Source of the data:", doi[tumor]), append = TRUE)

                } else {
                    createAlert(session, "tcgasearchmessage", "tcgasearchAlert", title = "Source of the data", style =  "success",
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
        updateCollapse(session, "collapseTCGA", open = "TCGA search results")
        output$tcgaSearchtbl <- renderDataTable({
            tumor <- isolate({input$tcgatumorClinicalFilter})
            type <- isolate({input$tcgaClinicalFilter})
            text.samples <- isolate({input$clinicalBarcode})
            tbl <- data.frame()
            withProgress( message = 'Search in progress',
                          detail = 'This may take a while...', value = 0, {
                              result = tryCatch({
                                  if(isolate({input$clinicalSearchType})){
                                      tbl <- rbind(tbl, TCGAquery_clinic(tumor = tumor, clinical_data_type = type))
                                  } else {
                                      samples <- parse.textarea.input(text.samples)
                                      tbl <- rbind(tbl, TCGAquery_clinic(samples = samples, clinical_data_type = type))
                                  }

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
                # Saving data
                if (isolate({input$saveClinicalRda}) || isolate({input$saveClinicalCsv})){
                    if(isolate({input$clinicalSearchType})){
                        filename <- paste0(tumor,"_clinic_",type,".rda")
                    } else {
                        filename <- paste0("samples_clinic_",type,gsub(" ","_",Sys.time()),".rda")
                    }
                    getPath <- parseDirPath(volumes, input$tcgafolder)
                    if (length(getPath) == 0) getPath <- "."
                    filename <- file.path(getPath,filename)
                    save.message <- ""
                    if (isolate({input$saveClinicalRda})){
                        save(tbl, file = filename)
                        save.message <-  paste0(save.message,"<br>File created: ",filename)
                    }
                    if (isolate({input$saveClinicalCsv})){
                        write.csv2(tbl, file = gsub("rda","csv",filename))
                        save.message <-  paste0(save.message,"<br>File created: ",gsub("rda","csv",filename))
                    }
                    createAlert(session, "tcgasearchmessage", "tcgasearchAlert", title = "File created", style =  "info",
                                content = save.message, append = TRUE)
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
    # Maf Search
    observeEvent(input$tcgaMafSearchBt, {

        output$tcgaSearchtbl <- renderDataTable({
            tumor <- isolate({input$tcgaMafTumorFilter})

            maf.files <- get.obj("maf.files")
            maf.files[,c(10,1,2,4,5,7)]
            tbl <- subset(maf.files,maf.files$Tumor == tumor)
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
            }
            return(tbl)
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

    observeEvent(input$mafDownloadBt , {
        maf.files <- get.obj("maf.files")
        # Dir to save the files
        getPath <- parseDirPath(volumes, input$tcgafolder)
        if (length(getPath) == 0) getPath <- "."

        withProgress(message = 'Download in progress',
                     detail = 'This may take a while...', value = 0, {
                         if(is.null(input$allRows)){
                             closeAlert(session, "tcgasearchAlert")
                             createAlert(session, "tcgasearchmessage", "tcgasearchAlert", title = "Error", style = "alert",
                                         content =  paste0("Please select which files will be downloaded"), append = TRUE)
                             req(input$allRows)
                         }
                         df <- data.frame(name = input$allRows[seq(5, length(input$allRows), 7)])
                         df <- maf.files[maf.files$Archive.Name %in% df$name,]

                         for (i in 1:nrow(df)) {
                             incProgress(1/nrow(df))
                             fout <- file.path(getPath,basename(df[1,]$Deploy.Location))
                             if (!file.exists(fout))  downloader::download(df[1,]$Deploy.Location,fout)
                         }})
        closeAlert(session, "tcgasearchAlert")
        createAlert(session, "tcgasearchmessage", "tcgasearchAlert", title = "Download completed", style = "success",
                    content =  paste0("Saved file: ",fout), append = FALSE)
    })
    #----------------------------------------------------------------------
    #                                         Summary plot
    #----------------------------------------------------------------------
    #-------------------------START controlling show/hide states -----------------

    observeEvent(input$summaryInputRb, {
        if(isolate({input$summaryInputRb}) == "platsample"){
            shinyjs::hide("tcgaSummaryExpFilter")
            shinyjs::show("tcgaSummarySamplestypeFilter")
            shinyjs::hide("summaryAddBarCount")

        } else {
            shinyjs::show("tcgaSummaryExpFilter")
            shinyjs::hide("tcgaSummarySamplestypeFilter")
            shinyjs::show("summaryAddBarCount")
        }
    })

    #-------------------------END controlling show/hide states -----------------
    observeEvent(input$tcgaSummaryBt , {
        output$summary.plot <- renderPlot({

            closeAlert(session, "tcgaSummaryAlert")

            tumor <- isolate({input$tcgaSummaryTumorFilter})
            platform <- isolate({input$tcgaSummaryExpFilter})
            level <- isolate({input$tcgaSummaryLevelFilter})

            if(is.null(tumor)){
                createAlert(session, "tcgaSummaryMessage", "tcgaSummaryAlert", title = "Data input error", style =  "danger",
                            content = "Please select a tumor type", append = FALSE)
                return(NULL)
            }

            if(isolate({input$summaryInputRb}) == "platsample"){


                result = tryCatch({
                    x <- TCGAquery(tumor = tumor,level=level)
                }, warning = function(w) {
                    next
                }, error = function(e) {
                    next
                })
                x <- TCGAquery(tumor = tumor,level=level)
                if(nrow(x)==0) next
                patient <- unique(substr(unlist(stringr::str_split(x$barcode,",")),1,15))

                samplesType <- isolate({input$tcgaSummarySamplestypeFilter})
                if (length(samplesType) > 0) {
                    patient <- unlist(lapply(samplesType,function(type){
                        patient[grep(type,substr(patient,14,15))]
                    }))
                }

                platform <- unique(x$Platform)
                df <- as.data.frame(matrix(0,nrow=length(patient),ncol=length(platform)+1))
                colnames(df) <- c("patient", platform)
                df$patient <- patient
                for (i in patient){
                    idx <- grep(i,x$barcode)
                    plat <- x[idx,"Platform"]
                    for (j in plat){
                        df[df$patient == i,j] <- 1
                    }
                }
                if(nrow(df) == 0) {
                    createAlert(session, "tcgaSummaryMessage", "tcgaSummaryAlert", title = "Results not found", style = "danger",
                                content =  paste0("No results found"), append = FALSE)
                    return(NULL)
                }
                df$Type <- tcga.code[substr(df$patient,14,15)]
                withProgress(message = 'Creating plot',
                             detail = 'This may take a while...', value = 0, {
                                 upset(df, nsets = length(platform),
                                       number.angles = 0,
                                       nintersects = 100,
                                       point.size = 3, name.size = 12,
                                       line.size = 1,
                                       mainbar.y.label = "Platform Intersections",
                                       sets.x.label = "Samples Per Platform",
                                       order.by = "freq",
                                       decreasing = T,
                                       #group.by = "sets",
                                       sets.bar.color = isolate({input$summarySetsBarColor}))
                             })
            } else {
                if(is.null(platform)){
                    createAlert(session, "tcgaSummaryMessage", "tcgaSummaryAlert", title = "Data input error", style =  "danger",
                                content = "Please select at least one platform", append = FALSE)
                    return(NULL)
                }
                not.found <- c()
                tbl <- data.frame()
                for(i in platform){
                    for(j in tumor){
                        x <- TCGAquery(tumor = j,platform = i,level=level)
                        if(!is.null(x)){
                            patient <- unique(substr(unlist(stringr::str_split(x$barcode,",")),1,15))
                            if(nchar(patient[1]) < 15) next
                            type <- tcga.code[substr(patient,14,15)]
                            tab <- table(substr(patient,14,15))
                            names(tab) <- tcga.code[names(tab)]
                            tab <- tab[!is.na(names(tab))]
                            tab <- (as.data.frame(tab))
                            tab$type <- rownames(tab)
                            tab$Platform <- i
                            tab$Tumor <- j
                            colnames(tab) <- c("Freq","type","Platform","Tumor")

                            if(nrow(tbl) == 0){
                                tbl <- tab
                            } else {
                                tbl <- rbind(tbl,tab)
                            }
                        } else {
                            not.found <- c(not.found, paste0("<li>Tumor: ",j,"  | Platform: ",i,"  | Level: ", level,"</li>"))
                        }
                    }
                }
                if(length(not.found) > 0) {
                    createAlert(session, "tcgaSummaryMessage", "tcgaSummaryAlert", title = "Results not found", style = "danger",
                                content =  paste0("These results were not found","<br><ul>", paste(not.found, collapse = ""),"</ul>"), append = FALSE)
                }
                if(nrow(tbl) == 0) {
                    return(NULL)
                }
                p <- ggplot(tbl, aes(x=type,y = Freq,fill=type)) + geom_bar(stat="identity") +
                    theme_bw() +theme(
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        axis.line.x = element_line(colour = "black"),
                        axis.line.y = element_line(colour = "black"),
                        legend.key = element_rect(colour = 'white'),
                        legend.justification=c(1,1),
                        legend.position=c(1,1),
                        text = element_text(size=16),
                        axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Type of sample") +
                    scale_fill_brewer(palette="Set1") + guides(fill=FALSE)
                #facet_wrap(~Platform + Tumor, ncol = isolate({input$summaryncol}))
                if(isolate({input$summaryAddBarCount})){
                    p <- p + geom_text(aes(vjust = 0, nudge_y = 0.7,label = Freq), size = 4)
                }
                p <- p + facet_grid(Platform ~ Tumor) +
                    theme(strip.background = element_rect(colour="navy", fill="navy",
                                                          size=1.5, linetype="solid"), strip.text=element_text(face= "bold",colour = "white"))
                plot(p)
            }
        })

    })
    observeEvent(input$tcgaSummaryBt , {
        updateCollapse(session, "collapseTCGAsummary", open = "Summary")
        output$tcgaSummary <- renderUI({
            plotOutput("summary.plot", width = paste0(isolate({input$summarywidth}), "%"), height = isolate({input$summaryheight}))
        })})



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
            shinyjs::show("oncoGenes")
        } else {
            shinyjs::show("oncoGenesTextArea")
            shinyjs::hide("oncoGenes")

        }
    })

    #-------------------------END controlling show/hide states -----------------

    shinyFileChoose(input, 'maffile', roots=volumes, session=session, restrictions=system.file(package='base'))
    shinyFileChoose(input, 'mafAnnotation', roots=volumes, session=session, restrictions=system.file(package='base'))
    annotation.maf <- function(){
        inFile <- input$mafAnnotation
        if (is.null(inFile)) return(NULL)
        file  <- as.character(parseFilePaths(volumes, input$mafAnnotation)$datapath)
        se <- get(load(file))

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
    observeEvent(input$oncoprintPlot , {

        output$oncoploting <- renderPlot({
            mut <- isolate({mut()})
            annotation <- isolate({annotation.maf()})
            cols <- isolate({input$mafAnnotationcols})
            rm.empty.cols <- isolate({input$oncoRmCols})
            show.col.names <- isolate({input$oncoShowColsNames})
            textarea <- isolate({input$oncoGenesTextArea})
            show.row.barplot <- isolate({input$oncoShowRowBarplot})

            if(isolate({input$oncoInputRb}) == "text" & !is.null(textarea)){
                genes <- toupper(parse.textarea.input(textarea))
                not.found <- genes[!(genes %in% mut$Hugo_Symbol)]
                if(length(not.found) > 0){
                    closeAlert("eaAlert")
                    createAlert(session, "oncomessage", "oncoAlert", title = "Data input error", style =  "danger",
                                content = paste0("Sorry, I cant't find these genes: ", not.found), append = FALSE)
                    genes <-  genes[genes %in% mut$Hugo_Symbol]
                }
            } else {
                genes <- isolate({input$oncoGenes})
            }

            if(is.null(genes)){
                createAlert(session, "oncomessage", "oncoAlert", title = "Error", style =  "danger",
                            content = "Please select the genes (max 50)", append = TRUE)
            } else if( length(genes) > 50){
                createAlert(session, "oncomessage", "oncoAlert", title = "Error", style =  "danger",
                            content = "The limit of the genes is 50", append = TRUE)
                return(NULL)
            }

            if(is.null(cols)) {
                annotation <- NULL
            } else if("bcr_patient_barcode" %in% colnames(annotation)) {
                annotation <- annotation[,c("bcr_patient_barcode",cols)]
            } else {
                createAlert(session, "oncomessage", "oncoAlert", title = "Error", style =  "danger",
                            content = "I couldn't find the bcr_patient_barcode column in the annotation", append = TRUE)
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

                             create.oncoprint(mut=mut,genes=genes,annotation = annotation,
                                              annotation.position=isolate(input$mafAnnotationpos),
                                              rm.empty.columns = rm.empty.cols,
                                              show.column.names = show.col.names,
                                              show.row.barplot = show.row.barplot,
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
        updateSelectizeInput(session, 'mafAnnotationcols', choices = as.character(colnames(annotation.maf())), server = TRUE)
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

    shinyFileChoose(input, 'volcanofile', roots=volumes, session=session, restrictions=system.file(package='base'))

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
                         df <- read.csv2(file,header = T)
                     })
        return(df)
    })

    observeEvent(input$volcanofile, {
        file  <- basename(as.character(parseFilePaths(volumes, input$volcanofile)$datapath))
        print(file)
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
    observeEvent(input$volcanoPlotBt , {
        output$volcano.plot <- renderPlot({

            # read csv file with results
            data <- isolate({volcanodata()})
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
                withProgress(message = 'Creating plot',
                             detail = 'This may take a while...', value = 0, {
                                 TCGAVisualize_volcano(x = data[,diffcol],
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
                                                       filename = NULL)
                             })
            } else {

                label <- c("Not Significant",
                           "Upregulated",
                           "Downregulated")
                label[2:3] <-  paste(label[2:3], "in", group2)
                if(isolate({input$volcanoNames})) names <- as.character(data$mRNA)

                withProgress(message = 'Creating plot',
                             detail = 'This may take a while...', value = 0, {
                                 TCGAVisualize_volcano(x = data$logFC,
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
                                                       filename = NULL)
                             })
            }
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
                                                   p.cut = isolate({input$dmrpvalue}),
                                                   diffmean.cut = isolate({input$dmrthrsld}),
                                                   cores = isolate({input$dmrcores}))
                             message <- paste0(message,"<li>DMR_results_", isolate({input$dmrgroupCol}), "_", group1, "_", group2, "_",
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
    shinyFileChoose(input, 'dmrfile', roots=volumes, session=session, restrictions=system.file(package='base'))
    shinyFileChoose(input, 'meanmetfile', roots=volumes, session=session, restrictions=system.file(package='base'))
    shinyFileChoose(input, 'heatmapfile', roots=volumes, session=session, restrictions=system.file(package='base'))
    shinyFileChoose(input, 'heatmapresultsfile', roots=volumes, session=session, restrictions=system.file(package='base'))

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
            if(!is.null(dmrdata())) as.character(colnames(colData(dmrdata())))
        }, server = TRUE)
    })
    observe({
        updateSelectizeInput(session, 'meanmetsubgroupCol', choices = {
            if(!is.null(meandata())) as.character(colnames(colData(meandata())))
        }, server = TRUE)
    })
    observe({
        updateSelectizeInput(session, 'meanmetgroupCol', choices = {
            if(!is.null(meandata())) as.character(colnames(colData(meandata())))
        }, server = TRUE)
    })




    observeEvent(input$meanmetPlot , {

        output$mean.plotting <- renderPlot({
            closeAlert(session, "meanmetAlert")
            jitter <- isolate({input$meanmetplotjitter})
            sort <- isolate({input$meanmetsort})
            angle <- isolate({input$meanmetAxisAngle})
            data <- meandata()

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

            withProgress(message = 'Creating plot',
                         detail = 'This may take a while...', value = 0, {
                             if(is.null(sort)){
                                 TCGAvisualize_meanMethylation(data=data,
                                                               groupCol=group,
                                                               subgroupCol=subgroup,
                                                               filename = NULL,
                                                               plot.jitter = jitter,
                                                               axis.text.x.angle = angle )
                             } else {
                                 TCGAvisualize_meanMethylation(data=data,
                                                               groupCol=group,
                                                               subgroupCol=subgroup,
                                                               filename = NULL,
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

    observeEvent(input$heatmapTypeInputRb, {
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
            if(input$heatmapProbesInputRb == "text"){
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
                         df <- read.csv2(file,header = T)
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
                            content = paste0("Please select the heatmapSortCol",append = FALSE))
                return(NULL)
            }

            if( isolate({input$heatmapTypeInputRb})=="met"){
                # ---------------- probes selection
                if(isolate({input$heatmapProbesInputRb}) == "Status"){
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
        } else if(input$tcgaEaInputRb == "Selection") {
            shinyjs::hide("eaGenesTextArea")
            shinyjs::show("eagenes")
        }
    })
    #----------------------- END controlling show/hide states -----------------
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
            } else {
                genes <- isolate({input$eagenes})
            }

            withProgress(message = 'Creating plot',
                         detail = 'This may take a while...', value = 0, {
                             ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor", genes)

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
                                                     nRGTab = genes,
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
        if(tools::file_ext(file)=="csv"){
            df <- read.csv2(file,header = T)
            rownames(df) <- df[,1]
            df[,1] <- NULL
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

            if(is.null(data)){
                closeAlert(session, "profileplotAlert")
                createAlert(session, "profileplotmessage", "profileplotAlert", title = "Missing data", style =  "danger",
                            content = paste0("Please select the data"), append = FALSE)
                return(NULL)
            }

            if(is.null(groupCol) || nchar(groupCol) == 0){
                closeAlert(session, "profileplotAlert")
                createAlert(session, "profileplotmessage", "profileplotAlert", title = "Missing group selection", style =  "danger",
                            content = paste0("Please select the group column"), append = FALSE)
                return(NULL)
            }

            if(is.null(subtypeCol) || nchar(subtypeCol) == 0){
                closeAlert(session, "profileplotAlert")
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
            closeAlert(session, "survivalAlert")
            createAlert(session, "survivalmessage", "survivalAlert", title = "Data input error", style =  "danger",
                        content = paste0("Sorry, but I'm expecting a csv or rda file, but I got a: ",
                                         tools::file_ext(file)), append = FALSE)
            return(NULL)
        }
        if(class(df)!= class(data.frame())){
            closeAlert(session, "survivalAlert")
            createAlert(session, "survivalmessage", "survivalAlert", title = "Data input error", style =  "danger",
                        content = paste0("Sorry, but I'm expecting a Data frame object, but I got a: ",
                                         class(df)), append = FALSE)
            return(NULL)
        }
        return(df)
    }

    observeEvent(input$survivalplotBt , {
        output$survival.plotting <- renderPlot({
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
                closeAlert(session, "survivalAlert")
                createAlert(session, "survivalmessage", "survivalAlert", title = "Missing data", style =  "danger",
                            content = paste0("Please select the data"), append = FALSE)
                return(NULL)
            }

            if(is.null(clusterCol) || nchar(clusterCol) == 0){
                closeAlert(session, "survivalAlert")
                createAlert(session, "survivalmessage", "survivalAlert", title = "Missing group", style =  "danger",
                            content = paste0("Please select group column"), append = FALSE)
                return(NULL)
            }

            if(length(unique(data[,clusterCol])) == 1){
                closeAlert(session, "survivalAlert")
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
    shinyFileChoose(input, 'survivalplotfile', roots=volumes, session=session, restrictions=system.file(package='base'))

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
                         exp[exp$logFC >= logFC.cut & exp$FDR <= fdr.cut,"status"] <- paste0("Upregulated in ", g2)
                         exp[exp$logFC <= -logFC.cut & exp$FDR <= fdr.cut,"status"] <- paste0("Downregulated in ", g2)
                     })

        out.filename <- paste0(paste("DEA_results",groupCol, g1, g2,"pcut",fdr.cut,"logFC.cut",logFC.cut,sep="_"),".csv")
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
            write.csv2(result$starburst, file = out.filename)
            createAlert(session, "starburstmessage", "starburstAlert", title = "Results saved", style =  "info",
                        content = paste0("Results saved in: ", out.filename), append = FALSE)
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
            se <- read.csv2(file,header = T)
            rownames(se) <- se[,1]
            se[,1] <- NULL
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
            se <- read.csv2(file,header = T)
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
    shinyFileChoose(input, 'starburstmetfile', roots=volumes, session=session, restrictions=system.file(package='base'))
    shinyFileChoose(input, 'starburstexpfile', roots=volumes, session=session, restrictions=system.file(package='base'))


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
    shinyFileChoose(input, 'elmermeefile', roots=volumes, session=session, restrictions=system.file(package='base'))
    shinyFileChoose(input, 'elmerresultsfile', roots=volumes, session=session, restrictions=system.file(package='base'))


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
                             Sig.probes <- get.diff.meth(mee,
                                                         cores=isolate({input$elmercores}),
                                                         dir.out =dir.out,
                                                         diff.dir=j,
                                                         pvalue = isolate({input$elmermetpvalue}))

                             #-------------------------------------------------------------
                             # Step 3.2: Identify significant probe-gene pairs            |
                             #-------------------------------------------------------------
                             # Collect nearby 20 genes for Sig.probes
                             nearGenes <- GetNearGenes(TRange=getProbeInfo(mee, probe=Sig.probes$probe),
                                                       cores=isolate({input$elmercores}),
                                                       geneAnnot=getGeneInfo(mee),
                                                       geneNum = isolate({input$elmergetpairNumGenes}))

                             Hypo.pair <- get.pair(mee=mee,
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
                             incProgress(1/2, detail = paste0('Analysis in direction completed: ',j))
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

}
