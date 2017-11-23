
getMafTumors <- function(){
    root <- "https://gdc-api.nci.nih.gov/data/"
    maf <- fread("https://gdc-docs.nci.nih.gov/Data/Release_Notes/Manifests/GDC_open_MAFs_manifest.txt",
                 data.table = FALSE, verbose = FALSE, showProgress = FALSE)
    tumor <- unlist(lapply(maf$filename, function(x){unlist(str_split(x,"\\."))[2]}))
    proj <- TCGAbiolinks:::getGDCprojects()

    disease <-  gsub("TCGA-","",proj$project_id)
    idx <- grep("disease_type",colnames(proj))
    names(disease) <-  paste0(proj[[idx]], " (",proj$project_id,")")
    disease <- sort(disease)
    ret <- disease[disease %in% tumor]
    return(ret)
}
observe({
    # Created with getMafTumors
    updateSelectizeInput(session, 'tcgaMafTumorFilter', choices =  maf.tumor, server = TRUE)
})

observeEvent(input$tcgaMafSearchBt, {
    output$tcgaMutationtbl <- DT::renderDataTable({
        tumor <- isolate({input$tcgaMafTumorFilter})
        getPath <- parseDirPath(get.volumes(isolate({input$workingDir})), input$workingDir)
        if (length(getPath) == 0) getPath <- paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI")

        withProgress(message = 'Download in progress',
                     detail = 'This may take a while...', value = 0, {
                         tbl <- GDCquery_Maf(tumor, directory = getPath,
                                             save.csv =  isolate({input$saveMafcsv}),
                                             pipelines = isolate({input$tcgaMafPipeline}))
                         incProgress(1, detail = "Download completed")
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
        if(grepl("varscan", isolate({input$tcgaMafPipeline}), ignore.case = TRUE)) {
            workflow.type <- "VarScan2 Variant Aggregation and Masking"
        } else if( isolate({input$tcgaMafPipeline}) == "muse") {
            workflow.type <- "MuSE Variant Aggregation and Masking"
        } else if( isolate({input$tcgaMafPipeline}) == "somaticsniper") {
            workflow.type <- "SomaticSniper Variant Aggregation and Masking"
        } else if(grepl("mutect", isolate({input$tcgaMafPipeline}), ignore.case = TRUE)) {
            workflow.type <-  "MuTect2 Variant Aggregation and Masking"
        } else {
            stop("Please select the pipeline argument (muse, varscan2, somaticsniper, mutect2)")
        }
        query <- GDCquery(paste0("TCGA-",tumor),
                          data.category = "Simple Nucleotide Variation",
                          data.type = "Masked Somatic Mutation",
                          workflow.type = workflow.type)

        fout <- file.path(getPath, getResults(query)$file_name)
        save(tbl,file = gsub("maf.gz","rda",fout))
        closeAlert(session, "tcgaMutationAlert")
        if(isolate({input$saveMafcsv})) {
            createAlert(session, "tcgaMutationmessage", "tcgaMutationAlert", title = "Download completed", style = "success",
                        content =  paste0("Saved file: <br><ul>",fout,"</ul><ul>", gsub(".gz",".csv",fout), "</ul>"), append = FALSE)
        } else {
            createAlert(session, "tcgaMutationmessage", "tcgaMutationAlert", title = "Download completed", style = "success",
                        content =  paste0("Saved file: ", gsub("maf.gz","rda",fout)), append = FALSE)
        }
        return(createTable(tbl))
    })
})


