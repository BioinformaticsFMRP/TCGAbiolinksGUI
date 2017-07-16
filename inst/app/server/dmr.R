##----------------------------------------------------------------------
#                             DMR analysis
##----------------------------------------------------------------------

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Handling visibility
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
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

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# UPDATING FIELDS AFTER DATA INPUT
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
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
            as.character(colnames(data))
        }
    }, server = TRUE)
})
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Analysis
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observeEvent(input$dmrAnalysis , {
    getPath <- parseDirPath(get.volumes(isolate({input$workingDir})), input$workingDir)
    if (length(getPath) == 0) getPath <- paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI")
    groups <- t(combn(isolate({input$dmrgroups}),2))
    print(groups)
    # read the data from the downloaded path
    # prepare it
    se <- isolate({dmrdata()})
    # Removes probes with all NA
    se <- subset(se,subset = (rowSums(!is.na(assay(se))) > 0))
    for(i in 1:nrow(groups)) {
        group1 <- groups[i,1]
        group2 <- groups[i,2]
        group1.col <- gsub("[[:punct:]]| ", ".", group1)
        group2.col <- gsub("[[:punct:]]| ", ".", group2)
        statuscol <- paste("status",group2.col,group1.col,sep = ".")
        results <- NULL
        withProgress(message = 'DMR analysis in progress',
                     detail = paste(group1," vs ", group2), value = 0, {
                         message <- "<br>Saving the results also in a csv file:<ul>"
                         if(!statuscol %in% colnames(values(se))){
                             step <- 1000
                             n <- nrow(se)
                             for(j in 0:floor(n/step)){
                                 end <- ifelse(((j + 1) * step) > n, n,((j + 1) * step))
                                 results <- rbind(results,
                                                  tryCatch({
                                                      values(TCGAanalyze_DMR(data = se[((j * step) + 1):end,],
                                                                             groupCol = isolate({input$dmrgroupCol}),
                                                                             group1 = group1,
                                                                             group2 = group2,
                                                                             plot.filename = FALSE,
                                                                             save = FALSE,
                                                                             calculate.pvalues.probes = isolate({input$dmrPvalues}),
                                                                             p.cut = isolate({input$dmrpvalue}),
                                                                             diffmean.cut = isolate({input$dmrthrsld}),
                                                                             cores = isolate({input$dmrcores})))
                                                  }, error = function(e) {return(NULL)
                                                  }))
                                 incProgress(1/(ceiling(n/step) + 1), detail = paste("Completed ", j + 1, " of ",ceiling(n/step)))
                             }
                             if(!is.null(results)){
                                 rownames(results) <- results[,grep("probeID|Composite.Element.REF",colnames(results))]
                                 results <- results[,!(colnames(results) %in% colnames(values(se)))]
                             }
                             if(isolate({input$dmrPvalues}) == "all"){
                                 values(se) <- cbind(values(se),results)
                             } else {
                                 se <- se[rownames(se) %in% rownames(results),]
                                 values(se) <- cbind(values(se),results)

                             }
                             incProgress(1/(ceiling(n/step) + 1), detail = "Saving results")
                         }
                         se <- TCGAanalyze_DMR(data = se,
                                               groupCol = isolate({input$dmrgroupCol}),
                                               group1 = group1,
                                               group2 = group2,
                                               save = TRUE,
                                               calculate.pvalues.probes = isolate({input$dmrPvalues}),
                                               save.directory = getPath,
                                               plot.filename = paste0("DMR_volcano_",group1,"_vs_",group2,".pdf"),
                                               p.cut = isolate({input$dmrpvalue}),
                                               diffmean.cut = isolate({input$dmrthrsld}),
                                               cores = isolate({input$dmrcores}))

                         message <- paste0(message,"<li>DMR_results_",
                                           file.path(getPath,
                                                     paste0(gsub("_",".",isolate({input$dmrgroupCol})),
                                                            "_", gsub("_",".",group1), "_", gsub("_",".",group2), "_",
                                                            "pcut_",isolate({input$dmrpvalue}), "_",
                                                            "meancut_",isolate({input$dmrthrsld}),".csv")),"</li>")
                         file  <- as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), input$dmrfile)$datapath)
                         if(!grepl("results",file)) file <- gsub(".rda","_results.rda",file)
                         save(se,file = file)
                         setProgress(1, detail = paste("Saving completed"))
                     })
    }
    createAlert(session, "dmrmessage", "dmrAlert", title = "DMR completed", style =  "success",
                content = paste0("Summarized Experiment object with results saved in: ", file, message,"<ul>"),
                append = FALSE)
})
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# File selection
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observe({
    shinyFileChoose(input, 'dmrfile', roots=get.volumes(input$workingDir), session=session,
                    restrictions=system.file(package='base'),filetypes=c('', 'rda'))
})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Data input
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
dmrdata <-  reactive({
    inFile <- input$dmrfile
    if (is.null(inFile)) return(NULL)
    file  <- as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), input$dmrfile)$datapath)

    withProgress(message = 'Loading data',
                 detail = 'This may take a while...', value = 0, {
                     result.file <- gsub(".rda","_results.rda",file)
                     if(file.exists(result.file)) {
                         se <- get(load(result.file))
                     } else {
                         se <- get(load(file))
                     }
                     incProgress(1, detail = "Completed")
                 })
    if(class(se)!= class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
        createAlert(session, "dmrmessage", "dmrAlert", title = "Data input error", style =  "danger",
                    content = paste0("Sorry, but I'm expecting a Summarized Experiment object, but I got a: ",
                                     class(se)), append = FALSE)
        return(NULL)
    }
    return(se)
})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Table
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=

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
