# -------------------------------------------------
# Network inference
# -------------------------------------------------
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# File selection
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observe({
    shinyFileChoose(input, 'networkfile', roots=get.volumes(input$workingDir),
                    session=session, restrictions=system.file(package='base'),
                    filetypes=c('rda', 'RData'))
})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Network inference
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observeEvent(input$networkBt , {
    print("Start inference")

    data <- networkdata()

    # Input should be expression data with genes in columns and samples in rows
    data <- t(data)

    method <- isolate({input$networkInferenceMethodRb})
    withProgress(message = 'Start inference',
                 detail = 'This may take a while...', value = 0, {
                     myadj <- TCGAanalyze_networkInference(data, optionMethod=method)
                 })

    # Saving output
    out.filename <- paste0(paste("network_results_method", method,sep="_"),".rda")
    getPath <- parseDirPath(get.volumes(isolate({input$workingDir})), input$workingDir)
    if (length(getPath) == 0) getPath <- paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI")

    save(myadj, file = out.filename)

    createAlert(session, "networkmessage", "networkAlert", title = "Analysis completed", style =  "success",
                content = paste0("File saved in: ",out.filename), append = FALSE)
})


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Data input
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
networkdata <-  reactive({
    inFile <- input$networkfile
    if (is.null(inFile)) return(NULL)
    file  <- as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), input$networkfile)$datapath)

    withProgress(message = 'Loading data',
                 detail = 'This may take a while...', value = 0, {
                     data <- get(load(file))
                     incProgress(1, detail = "Completed")
                 })
    if(class(data) == class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
        data <- assays(data)[[2]]
    }

    # Check data input, should it be a data frame/matrix or summarizedExperiment
    if(class(data) != class(as(SummarizedExperiment(),"RangedSummarizedExperiment")) & class(data) != class(data.frame())){
        createAlert(session, "networkmessage", "networkAlert", title = "Data input error", style =  "danger",
                    content = paste0("Sorry, but I'm expecting either a Summarized Experiment object or a data frame/matrix, but I got a: ",
                                     class(data)), append = FALSE)
        return(NULL)
    }
    return(data)
})
