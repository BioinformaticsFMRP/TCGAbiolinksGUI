#-------------------------------------------------------------------------
#                            Summarized experiment edition
#-------------------------------------------------------------------------
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# DATA INPUT
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
seeditdata <- reactive({
    inFile <- input$seeditfile
    if(class(inFile) != "list") return(NULL)
    file  <- as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), input$seeditfile)$datapath)
    if(tools::file_ext(file)=="rda"){
        se <- get(load(file))
        return(se)
    }
    return(se)
})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# File selection
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observe({
    shinyFileChoose(input, 'seeditfile',
                    roots = get.volumes(input$workingDir),
                    session = session,
                    restrictions = system.file(package='base'),
                    filetypes = c('','rda'))
})

output$downloadData <- downloadHandler(
    filename = function() {
        file  <- basename(as.character(parseFilePaths(
            get.volumes(isolate({input$workingDir})), input$seeditfile)$datapath))
        paste0("samples-matrix-",gsub("rda",'csv',file))
    },
    content = function(con) {
        if(!is.null(seeditdata())) {
            data <- as.data.frame(colData(seeditdata()))
            file  <- as.character(parseFilePaths(
                get.volumes(isolate({input$workingDir})), input$seeditfile)$datapath)
            row.names <- rownames(data)
            write_csv(cbind(row.names,data), con)
        }
    }
)
observeEvent(input$seeditSEbt, {
    data <- seeditdata()
    closeAlert(session,"seeditAlert")
    inFile <- isolate({input$seuploadfile})
    if(class(inFile) != "list") return(NULL)
    sample.data <-   read_csv(inFile$datapath,col_names = TRUE)
    sample.data <- DataFrame(sample.data[match(sample.data$row.names,colnames(data)),])
    rownames(sample.data) <- sample.data$row.names
    sample.data$row.names <- NULL
    colData(data) <- sample.data
    getPath <- parseDirPath(get.volumes(isolate({input$workingDir})), input$workingDir)
    if (length(getPath) == 0) getPath <- paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI")
    filename <- file.path(getPath,isolate({input$seeditfilename}))
    withProgress(message = 'Saving new summarized Experiment object',
                 detail = 'This may take a while...', value = 0, {
                     save(data, file = filename)
                     incProgress(1, detail = "Completed")

                 })
    createAlert(session,
                "seeditAlert",
                "seeditAlertMessage",
                title = "File created",
                style =  "success",
                content = paste0("Saved as: ",filename), append = FALSE)
})
observeEvent(input$seeditfile, {
    updateTextInput(session,
                    "seeditfilename",
                    value = paste0("edited_",
                                   basename(as.character(parseFilePaths(
                                       get.volumes(isolate({input$workingDir})),
                                       input$seeditfile)$datapath))))

    updateCollapse(session, "collapseseeEdit", open = "Sample data matrix")

    output$seeditassaytb <- DT::renderDataTable({
        x <- as.data.frame(assay(seeditdata()))
        row.names <- rownames(x)
        x <- cbind(row.names,x)
        return(createTable(x))
    })

    output$seeditfeaturetb <- DT::renderDataTable({
        return(createTable(as.data.frame(rowRanges(seeditdata()))))
    })

    output$seeditsampletb <- DT::renderDataTable({
        return(createTable(as.data.frame(colData(seeditdata()))))
    })
})
