#------------------------------------------------
# Starburst
# -----------------------------------------------

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Handling visibility
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observeEvent(input$starburstNames, {
    toggle("starburstNamesFill")
})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# UPDATING FIELDS AFTER DATA INPUT
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# get DMR result name and update the mean cut and pcut
observeEvent(input$starburstmetfile, {
    file  <- basename(as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), input$starburstmetfile)$datapath))
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
    file  <- basename(as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), input$starburstexpfile)$datapath))
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

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Plot
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=

# Main function
# Main function
starburst <-  reactive({
    trigger <- input$starburstPlot
    closeAlert(session,"starburstAlert")
    exp <- isolate({result.dea.data()})
    met <- isolate({result.dmr.data()})
    if(is.null(exp)){
        createAlert(session, "starburstmessage", "starburstAlert", title = "Missing data", style =  "danger",
                    content = paste0("Please select the differential expression results"), append = FALSE)
        return(NULL)
    }
    if(is.null(met)){
        createAlert(session, "starburstmessage", "starburstAlert", title = "Missing data", style =  "danger",
                    content = paste0("Please select the differential DNA methylation results"), append = FALSE)
        return(NULL)
    }
    logFC.cut <- isolate({input$starburstexpFC})
    exp.p.cut <- isolate({input$starburstexFDR})
    diffmean.cut <- isolate({input$starburstmetdiff})
    met.p.cut <- isolate({input$starburstmetFDR})

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


    file  <- basename(as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), isolate({input$starburstmetfile}))$datapath))
    if(length(file) > 0){
        file <- unlist(str_split(file,"_"))
        group1 <- file[4]
        group2 <- file[5]
    }
    file  <- basename(as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), isolate({input$starburstexpfile}))$datapath))
    if(length(file) > 0){
        file <- unlist(str_split(file,"_"))
        exp.group1 <- file[4]
        exp.group2 <- file[5]
    }
    # As the DNA methylation group has both relations A vs B and B vs A we will just change it to
    # be in the same relation of the DEA results. A vs B. the logFC is log2(B/A)
    if(gsub("[[:punct:]]| ", ".", group1) == gsub("[[:punct:]]| ", ".", exp.group2)
       & gsub("[[:punct:]]| ", ".", group2) == gsub("[[:punct:]]| ", ".", exp.group1)){
        aux <- group1
        group1 <- group2
        group2 <- aux
    }
    if(gsub("[[:punct:]]| ", ".", group1) == gsub("[[:punct:]]| ", ".", exp.group1)
       & gsub("[[:punct:]]| ", ".", group2) == gsub("[[:punct:]]| ", ".", exp.group2)){
        if(!any(grepl("ENSG",exp[1,]))) rownames(exp) <- exp$Gene_symbol
        result <- TCGAvisualize_starburst(met = met,
                                          exp = exp,
                                          group1 = group1,
                                          group2 = group2,
                                          color = colors,
                                          genome = isolate({input$starburstGenome}),
                                          met.platform = isolate({input$starburstMetPlatform}),
                                          names = names,
                                          names.fill = names.fill,
                                          exp.p.cut = exp.p.cut,
                                          met.p.cut = met.p.cut,
                                          diffmean.cut = diffmean.cut,
                                          logFC.cut = logFC.cut,
                                          return.plot = TRUE)
        out.filename <- paste0(paste("Starburst_results", group1, group2,
                                     "exp.p.cut", exp.p.cut, "logFC.cut", logFC.cut,
                                     "met.diffmean", diffmean.cut, "met.p.cut", met.p.cut,
                                     sep = "_"),".csv")
        if(isolate({input$starburstSave})) {
            write_csv(result$starburst, path = out.filename)
            createAlert(session, "starburstmessage", "starburstAlert", title = "Results saved", style =  "success",
                        content = paste0("Results saved in: ", out.filename), append = FALSE)
        }
        return(result)
    } else {
        createAlert(session, "starburstmessage", "starburstAlert", title = "Error", style =  "error",
                    content = paste0("The sames groups were not found in the results files"), append = FALSE)
    }
})
# -------------- Starburst plot
observeEvent(input$starburstPlot , {
    # validate input
    output$starburst.plot <- renderPlot({
        withProgress(message = 'Creating plot',
                     detail = 'This may take a while...', value = 0, {
                         aux <- starburst()

                         if(!is.null(aux)) {
                             incProgress(1, detail = "Rendering plot")
                             return(aux$plot)
                         }
                     })
    })})

observeEvent(input$starburstPlot , {
    updateCollapse(session, "collapsedea", open = "dea plots")
    output$starburstPlot <- renderUI({
        plotOutput("starburst.plot", width = paste0(isolate({input$starburstwidth}), "%"), height = isolate({input$starburstheight}))
    })})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# DATA INPUT
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
result.dea.data <-  reactive({
    inFile <- input$starburstexpfile
    if (is.null(inFile)) return(NULL)
    file  <- as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), input$starburstexpfile)$datapath)
    if(tools::file_ext(file)=="csv"){
        se <- as.data.frame(read_csv(file)); se$X1 <- NULL
    } else if(tools::file_ext(file)=="rda"){
        se <- get(load(file))
    }
    if(all(grepl("ENS", se$Gene_symbol))) {
        fromType <- "ENSEMBL"
        # In case we have ENSG
        eg = as.data.frame(bitr(se$Gene_symbol,
                                fromType=fromType,
                                toType="SYMBOL",
                                OrgDb="org.Hs.eg.db"))
        eg <- eg[!duplicated(eg[,fromType]),]
        colnames(se)[grep("Gene_symbol",colnames(se))] <- "Gene"
        colnames(eg) <- c("Gene","Gene_symbol")
        se <- merge(se,eg,by = "Gene")
    }

    if(class(se)!= class(data.frame())){
        createAlert(session, "deamessage", "deaAlert", title = "Data input error", style =  "danger",
                    content = paste0("Sorry, but I'm expecting a Data frame object, but I got a: ",
                                     class(se)), append = FALSE)
        return(NULL)
    }
    return(se)
})
result.dmr.data <-  reactive({
    inFile <- input$starburstmetfile
    if (is.null(inFile)) return(NULL)
    file  <- as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), input$starburstmetfile)$datapath)
    if(tools::file_ext(file)=="csv"){
        se <- as.data.frame(read_csv(file)); se$X1 <- NULL
    } else if(tools::file_ext(file)=="rda"){
        se <- get(load(file))
    }
    return(se)
})
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# File selection
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observe({
    shinyFileChoose(input, 'starburstmetfile', roots=get.volumes(input$workingDir), session=session,
                    restrictions=system.file(package='base'), filetypes=c('excel', 'csv'))
    shinyFileChoose(input, 'starburstexpfile', roots=get.volumes(input$workingDir), session=session,
                    restrictions=system.file(package='base'), filetypes=c('excel', 'csv'))
})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Table output with results
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=

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


output$savestarburstpicture <- downloadHandler(
    filename = function(){input$starburstPlot.filename},
    content = function(file) {
        if(tools::file_ext(input$starburstPlot.filename) == "png") {
            device <- function(..., width, height) {
                grDevices::png(..., width = 10, height = 10,
                               res = 300, units = "in")
            }
        } else if(tools::file_ext(input$starburstPlot.filename) == "pdf") {
            device <- function(..., width, height) {
                grDevices::pdf(..., width = 10, height = 10)
            }
        } else if(tools::file_ext(input$starburstPlot.filename) == "svg") {
            device <- function(..., width, height) {
                grDevices::svg(..., width = 10, height = 10)
            }
        }

        ggsave(file, plot =  starburst()$plot, device = device)
    })
