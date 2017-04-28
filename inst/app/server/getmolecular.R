#-------------------------------------------------------------------------
#                            TCGA Search
#-------------------------------------------------------------------------
#--------------------- START controlling show/hide states ----------------
observeEvent(input$clinicalIndexed, {
    if(input$clinicalIndexed){
        shinyjs::hide("tcgaClinicalFilter")
    } else {
        shinyjs::show("tcgaClinicalFilter")
    }
})
query.result <-  reactive({
    input$tcgaSearchBt # the trigger
    #------------------- STEP 1: Argument-------------------------
    # Set arguments for GDCquery, if value is empty we will set it to FALSE (same as empty)
    tumor <- isolate({input$tcgaProjectFilter})
    if(str_length(tumor) == 0)  tumor <- NULL

    data.category <- isolate({input$tcgaDataCategoryFilter})
    if(str_length(data.category) == 0)  data.category <- NULL

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

    # Close notifications
    if (!is.null(id)) {
        removeNotification(id)
        id <<- NULL
    }

    if(is.null(tumor)){
        id <<- showNotification("Input error: Please select a project", type =  "error", duration = NULL)
        return(NULL)
    }
    if(is.null(data.category)){
        id <<- showNotification("Input error: Please select a data.category", type =  "error", duration = NULL)
        return(NULL)
    }

    withProgress(message = 'Accessing GDC',
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
                         incProgress(1, detail = "Completed")
                         return(query)
                     }, error = function(e) {
                         id <<- showNotification("Search error: No results for this query", type =  "error", duration = NULL)
                         return(NULL)
                     })
                 })

    not.found <- c()
    tbl <- data.frame()
    results <- query$results[[1]]


    #------------------- STEP 2: Clinical features-------------------------
    # In this step we will download indexed clinical data
    # Get the barcodes that repects user inputs
    # And select only the results that respects it

    clinical <- getClinical.info()

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
    if(any(!is.na(gender.idx),!is.na(race.idx),!is.na(vital.status.idx),!is.na(stage.idx))){
        idx <- apply(data.frame(gender.idx,race.idx,vital.status.idx,stage.idx),1,function(x)all(x,na.rm = TRUE))
        clinical <- clinical[idx,]
    }
    cut <- ifelse(grepl("TARGET",query$project),16,12)
    results <- results[substr(results$cases,1,cut) %in% clinical$bcr_patient_barcode,]
    query$results[[1]] <- results
    list(query,clinical)
})
#----------------------- END controlling show/hide states -----------------
observeEvent(input$tcgaSearchBt, {
    updateTextInput(session, "tcgafilename", label =  "File name",
                    value = paste0(
                        paste(isolate({input$tcgaProjectFilter}),
                              gsub(" ","_",isolate({input$tcgaDataCategoryFilter})),
                              gsub(" ","_",isolate({input$tcgaDataTypeFilter})),
                              ifelse(isolate({input$tcgaDatabase}),"hg19","hg38"),
                              sep = "_"),".rda"))
    updateCollapse(session, "collapseTCGA", open = "GDC search results")

    output$queryresutlstable <- renderDataTable({
        results <- getResults(query.result()[[1]])
        if(!is.null(results)) return(as.data.frame(results))
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

    output$tcgaview <- renderGvis({

        results <- getResults(query.result()[[1]])
        clinical <- query.result()[[2]]
        if(any(duplicated(results$cases))) id <<- showNotification("Search: There are more than one file for the same case.", type =  "warning", duration = NULL)

        #------------------- STEP 3: Plot-------------------------
        # Plot informations to help the user visualize the data that will
        # be downloaded
        if(is.null(results)) return(NULL)
        legacy <- isolate({as.logical(input$tcgaDatabase)})

        data.type <- gvisPieChart(as.data.frame(table(results$data_type)),
                                  options=list( title="Data type"))
        tissue.definition <- gvisPieChart(as.data.frame(table(results$tissue.definition)),
                                          options=list( title="Tissue definition"))
        experimental_strategy  <- gvisPieChart(as.data.frame(table(results$experimental_strategy )),
                                               options=list( title="Experimental strategy"))
        analysis.workflow_type <- gvisPieChart(as.data.frame(table(paste0(
            results$analysis_workflow_type,"(", results$analysis_workflow_link, ")"))),
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
    # Close notifications
    if (!is.null(id)) {
        removeNotification(id)
        id <<- NULL
    }

    query <- query.result()[[1]]
    results <- query$results[[1]]
    # Dir to save the files
    getPath <- parseDirPath(get.volumes(isolate({input$workingDir})), input$workingDir)
    if (length(getPath) == 0) getPath <- paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI")
    filename <- file.path(getPath,isolate({input$tcgafilename}))
    withProgress(message = 'Download in progress',
                 detail = 'This may take a while...', value = 0, {
                     trash = tryCatch({
                         n <- nrow(results)
                         step <- 20
                         for(i in 0:ceiling(n/step - 1)){
                             query.aux <- query
                             end <- ifelse(((i + 1) * step) > n, n,((i + 1) * step))
                             query.aux$results[[1]] <- query.aux$results[[1]][((i * step) + 1):end,]
                             GDCdownload(query.aux, method = "api",directory = getPath)
                             incProgress(1/ceiling(n/step), detail = paste("Completed ", i + 1, " of ",ceiling(n/step)))
                         }
                     }, error = function(e) {
                         id <<- showNotification("Downloead error: Error while downloading the files. Please try again", type =  "warning", duration = NULL)
                         return(NULL)
                     })
                 })
    withProgress(message = 'Prepare progress',
                 detail = 'This may take a while...', value = 0, {
                     trash = tryCatch({
                         genes <- NULL
                         if(isolate({input$addGistic})) genes <- isolate({input$gisticGenes})
                         trash <- GDCprepare(query, save = TRUE,
                                             save.filename = filename,
                                             summarizedExperiment = as.logical(isolate({input$prepareRb})),
                                             directory = getPath,
                                             mut.pipeline = isolate({input$tcgaPrepareMutPipeline}),
                                             add.gistic2.mut = genes)
                         if(!is.null(genes)) {
                             aux <- gsub(".rda","_samples_information.csv",filename)
                             write_csv(x = as.data.frame(colData(trash)),path  = aux )
                             filename <- c(filename, aux)
                         }
                     }, error = function(e) {
                         id <<- showNotification("Error while preparing the files", type =  "error", duration = NULL)
                         print(e)
                         return(NULL)
                     })
                     createAlert(session, "tcgasearchmessage", "tcgaAlert", title = "Prepare completed", style =  "success",
                                 content =  paste0("Saved in: ", "<br><ul>", paste(filename, collapse = "</ul><ul>"),"</ul>"), append = FALSE)
                 })

})

getClinical.info <-  reactive({
    project <- input$tcgaProjectFilter
    if(is.null(project) || str_length(project) == 0) return(NULL)

    baseURL <- "https://gdc-api.nci.nih.gov/cases/?"
    options.pretty <- "pretty=true"
    options.expand <- "expand=diagnoses,demographic"
    option.size <- paste0("size=", TCGAbiolinks:::getNbCases(project, "Clinical"))
    files.data_category <- "Clinical"
    options.filter <- paste0("filters=", URLencode("{\"op\":\"and\",\"content\":[{\"op\":\"in\",\"content\":{\"field\":\"cases.project.project_id\",\"value\":[\""),
                             project, URLencode("\"]}},{\"op\":\"in\",\"content\":{\"field\":\"files.data_category\",\"value\":[\""),
                             files.data_category, URLencode("\"]}}]}"))

    withProgress(message = 'Loading clinical data',
                 detail = 'This may take a while...', value = 0, {
                     json <- jsonlite::fromJSON(paste0(baseURL, paste(options.pretty, options.expand,
                                                                      option.size, options.filter, sep = "&")), simplifyDataFrame = TRUE)
                 })
    results <- json$data$hits
    diagnoses <- rbindlist(results$diagnoses, fill = TRUE)
    diagnoses$submitter_id <- results$submitter_id
    results$demographic$submitter_id <- results$submitter_id
    df <- merge(diagnoses, results$demographic, by = "submitter_id", all = TRUE)
    df$bcr_patient_barcode <- df$submitter_id
    df$disease <- gsub("TCGA-|TARGET-", "", project)
    setDF(df)
    return(df)
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
    updateSelectizeInput(session, 'tcgaProjectFilter', choices =  getTCGAdisease(), server = TRUE)
    updateSelectizeInput(session, 'tcgatumorClinicalFilter', choices =  getTCGAdisease(), server = TRUE)
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
    input$tcgaProjectFilter
    tryCatch({
        clin <- getClinical.info()
        if(!is.null(clin)){
            updateSelectizeInput(session, 'tcgaClinicalGenderFilter', choices =  unique(clin$gender), server = TRUE)
            updateSelectizeInput(session, 'tcgaClinicalVitalStatusFilter', choices =  unique(clin$vital_status), server = TRUE)
            updateSelectizeInput(session, 'tcgaClinicalRaceFilter', choices =  unique(clin$race), server = TRUE)
            updateSelectizeInput(session, 'tcgaClinicalTumorStageFilter', choices =  unique(clin$tumor_stage), server = TRUE)
        }
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


observeEvent(input$prepareRb, {
    if(input$prepareRb) {
        shinyjs::show("addGistic")
        if(isolate({input$addGistic})){
            shinyjs::show("gisticGenes")
            shinyjs::show("tcgaPrepareMutPipeline")
        } else {
            shinyjs::hide("gisticGenes")
            shinyjs::hide("tcgaPrepareMutPipeline")
        }
    } else {
        shinyjs::hide("addGistic")
        shinyjs::hide("gisticGenes")
        shinyjs::hide("tcgaPrepareMutPipeline")
    }
})

observeEvent(input$addGistic, {
    if(input$addGistic) {
        shinyjs::show("gisticGenes")
    } else {
        shinyjs::hide("gisticGenes")
    }
})

observe({
    updateSelectizeInput(session, 'tcgaDataCategoryFilter', choices =  getDataCategory(as.logical(input$tcgaDatabase)), server = TRUE)
    updateSelectizeInput(session, 'gisticGenes', choices = as.character(sort(TCGAbiolinks:::gene.list)), server = TRUE)
})
