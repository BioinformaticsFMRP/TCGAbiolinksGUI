
observeEvent(input$tcgaClinicalBt, {
    updateCollapse(session, "collapseTCGAClinical", open = "Clinical data", close = "Description of the data")
    output$tcgaClinicaltbl <- DT::renderDataTable({
        closeAlert(session, "tcgaClinicalAlert")
        project <- isolate({input$tcgatumorClinicalFilter})
        type <- isolate({input$tcgaClinicalTypeFilter})
        parser <- isolate({input$tcgaClinicalFilter})
        text.samples <- isolate({input$clinicalBarcode})
        tbl <- data.frame()

        getPath <- parseDirPath(get.volumes(isolate({input$workingDir})), input$workingDir)
        if (length(getPath) == 0) getPath <- paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI")

        withProgress( message = 'Download in progress',
                      detail = 'This may take a while...', value = 0, {
                          result = tryCatch({
                              if(isolate({input$clinicalIndexed})){
                                  tbl <- rbind(tbl, GDCquery_clinic(project = project, type = type))
                              } else {
                                  if(grepl("clinical", type, ignore.case = TRUE)) type <- "Clinical"
                                  if(grepl("biospecimen", type, ignore.case = TRUE)) type <- "Biospecimen"
                                  query <- GDCquery(project = project, data.category = type)

                                  result = tryCatch({
                                      GDCdownload(query, directory = getPath)
                                  } , error = function(e) {
                                      createAlert(session, "tcgaClinicalmessage", "tcgaClinicalAlert", title = "No results found", style =  "error",
                                                  content = "There was a problem to download the data. Please try again later.", append = FALSE)

                                  })
                                  incProgress(1/2, message = "Reading XML files",detail = paste(" Parsing"))
                                  clinical <- GDCprepare_clinic(query, parser, directory = getPath)
                                  tbl <- rbind(tbl, clinical)
                              }
                              incProgress(1, detail = "Completed")
                              return(tbl)
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
                    if(grepl("clinical",type, ignore.case = TRUE))   tbl$treatments <- NULL
                    write_csv(tbl,path = gsub("rda","csv",filename))
                    save.message <-  paste0(save.message,"File created: ",gsub("rda","csv",filename))
                }
                createAlert(session, "tcgaClinicalmessage", "tcgaClinicalAlert", title = "File created", style =  "info",
                            content = save.message, append = TRUE)
            }
            return(createTable(tbl))
        }

    })
})


observeEvent(input$tcgaClinicalTypeFilter, {
    if(isolate({input$tcgaClinicalTypeFilter}) == "biospecimen") {
        options <- sort(c("protocol","admin","aliquot","analyte","bio_patient","sample", "portion", "slide"))
    } else {
        options <- sort(c("drug", "admin", "follow_up", "radiation", "patient", "stage_event", "new_tumor_event"))
    }
    updateSelectizeInput(session, 'tcgaClinicalFilter', choices = as.character(options), server = TRUE)
})
