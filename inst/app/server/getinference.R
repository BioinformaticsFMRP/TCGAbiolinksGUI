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
    output$tcgaNetworktbl <- renderDataTable({
        # Close in case it is open
        closeAlert(session, "networkAlert")
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

        vals <- unique(sort(myadj[upper.tri(myadj)],decreasing=TRUE))[1:10]

        tbl <- NULL
        for(i in 1:10){
            ind <-  which(myadj == vals[i], arr.ind = TRUE)[,2]
            for(j in 1:(length(ind)/2)) {
                aux <- data.frame("gene 1" = rownames(myadj)[ind[j]],  "gene 2" = rownames(myadj)[ind[j + 1]], "value" = vals[i])
                tbl <- rbind(tbl, aux)
            }
        }


        return(tbl)
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
    )
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

    # Check data input, should it be a data frame/matrix or summarizedExperiment
    if(class(data) != class(as(SummarizedExperiment(),"RangedSummarizedExperiment")) & class(data) != class(data.frame())){
        createAlert(session, "networkmessage", "networkAlert", title = "Data input error", style =  "danger",
                    content = paste0("Sorry, but I'm expecting either a Summarized Experiment object or a data frame/matrix, but I got a: ",
                                     class(data)), append = FALSE)
        return(NULL)
    }

    if(class(data) == class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
        data <- assay(data)
    }

    return(data)
})
