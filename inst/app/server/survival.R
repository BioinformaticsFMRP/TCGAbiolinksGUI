#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Survival plot
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Handling visibility
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# For survival we have the analysis by group or by the gene expression level
# We will hide/show the fields based on the type of analysis
observeEvent(input$survivalbyGene , {
    if(input$survivalbyGene){
        shinyjs::show("survivalPercent")
        shinyjs::show("survivalGene")
        shinyjs::hide("survivalplotgroup")
    } else {
        shinyjs::hide("survivalPercent")
        shinyjs::hide("survivalGene")
        shinyjs::show("survivalplotgroup")
    }
})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# DATA INPUT
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
survivalplotdata <- function(){
    closeAlert(session, "survivalAlert")
    inFile <- input$survivalplotfile
    if (is.null(inFile)) return(NULL)
    file  <- as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), inFile)$datapath)
    if(tools::file_ext(file)=="csv"){
        df <- as.data.frame(read_csv(file, col_names = TRUE))
        if(ncol(df) == 1) df <- as.data.frame(read_csv(file, col_names = TRUE))
        rownames(df) <- df[,1]
        df[,1] <- NULL
    } else if(tools::file_ext(file)=="rda"){
        df <- get(load(file))
    } else {
        createAlert(session, "survivalmessage", "survivalAlert",
                    title = "Data input error", style =  "danger",
                    content = paste0("Sorry, but I'm expecting a csv or rda file, but I got a: ",
                                     tools::file_ext(file)), append = FALSE)
        return(NULL)
    }

    # Check required columns
    cols <- c("days_to_death","days_to_last_follow_up","vital_status")
    if(class(df) ==  class(data.frame())){
        if(!(all(cols %in% colnames(df)))){
            createAlert(session, "survivalmessage", "survivalAlert",
                        title = "Data input error", style =  "danger",
                        content = paste0("Sorry, but I'm expecting columns: ",
                                         paste(cols,collapse = ", ")), append = FALSE)
            return(NULL)
        }
    } else {
        if (!(all(cols %in% colnames(colData(df))))) {
            createAlert(session, "survivalmessage", "survivalAlert",
                        title = "Data input error", style =  "danger",
                        content = paste0("Sorry, but I'm expecting columns: ",
                                         paste(cols,collapse = ", ")), append = FALSE)
            return(NULL)
        }
    }
    return(df)
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# UPDATING FIELDS AFTER DATA INPUT
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Update groups to be used in the analysis
observe({
    data <- survivalplotdata()
    updateSelectizeInput(session, 'survivalplotgroup', choices = {
        if (!is.null(data)) {
            if (class(data) ==  class(data.frame())) choices <- as.character(colnames(data))
            if (class(data) !=  class(data.frame())) choices <- as.character(colnames(colData(data)))
            choices
        }
    }, server = TRUE)
})

# Update genes list to be used in the analysis of survival
observe({
    data <- survivalplotdata()
    if(class(data) == class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
        shinyjs::show("survivalbyGene")

        updateSelectizeInput(session, 'survivalGene', choices = {
            if(!is.null(data)) rownames(data)
        }, server = TRUE)
    } else {
        updateCheckboxInput(session, "survivalbyGene",value = FALSE)
        shinyjs::hide("survivalbyGene")

    }
})

observe({
    data <- survivalplotdata()
    updateSelectizeInput(session, 'survivalplotsubtype', choices = {
        if(!is.null(data)) as.character(colnames(data))
    }, server = TRUE)
})
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observeEvent(input$survivalplotgroup , {
    if(isolate({input$survivalplotgroup}) != "") {
        updateTextInput(session, "survivalplotLegend", value = isolate({input$survivalplotgroup}))
    }
})
observeEvent(input$survivalGene , {
    if(isolate({input$survivalGene}) != "") {
        updateTextInput(session, "survivalplotLegend", value = paste0("Expression of ",isolate({input$survivalGene})))
    }
})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Plot
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
survival.plot <-  reactive({
    input$survivalplotBt
    data <- isolate({survivalplotdata()})
    clusterCol <-  isolate({input$survivalplotgroup})
    legend <- isolate({input$survivalplotLegend})
    main <- isolate({input$survivalplotMain})
    xlim <- isolate({input$survivalplotLimit})
    if(xlim == 0) xlim <- NULL
    pvalue <- isolate({input$survivalplotPvalue})
    risk.table <- isolate({input$survivalplotRiskTable})
    conf.int <- isolate({input$survivalplotConfInt})

    # if Summarized Experiment
    if(isolate(input$survivalbyGene)){
        aux <- assay(data)[rownames(data) == isolate({input$survivalGene}),]
        colData(data)$level.group <- "Mid expression"
        min.cut <- max(sort(aux)[1:(length(aux) * isolate({input$survivalPercent}))])
        high.cut <- min(sort(aux, decreasing = T)[1:(length(aux) * isolate({input$survivalPercent}))])
        colData(data)[aux <= min.cut,"level.group"] <- "Low expression"
        colData(data)[aux >= high.cut,"level.group"] <- "High expression"
        clusterCol <- "level.group"
    }
    if(class(data) !=  class(data.frame())){
        data <- colData(data)
    }
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
                    content = paste0("Sorry, but I'm expecting at least two groups<br>",
                                     "Only this group found: ", unique(data[,clusterCol])), append = FALSE)
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
                                          xlim = xlim,
                                          pvalue = pvalue,
                                          risk.table=risk.table,
                                          conf.int=conf.int)

                 })
})

observeEvent(input$survivalplotBt , {
    output$survival.plotting <- renderPlot({
        closeAlert(session, "survivalAlert")
        survival.plot()
    })
})
observeEvent(input$survivalplotBt , {
    updateCollapse(session, "collapsesurvivalplot", open = "Survival plot")
    output$survivalplot <- renderUI({
        plotOutput("survival.plotting",
                   width = paste0(isolate({input$survivalwidth}), "%"),
                   height = isolate({input$survivalheight}))
    })})
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# File selection
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observe({
    shinyFileChoose(input, 'survivalplotfile',
                    roots=get.volumes(input$workingDir),
                    session=session,
                    restrictions=system.file(package='base'),
                    filetypes=c('', 'csv','rda'))
})

output$savesurvivalpicture <- downloadHandler(
    filename = function(){input$survivalPlot.filename},
    content = function(file) {
        if(tools::file_ext(input$survivalPlot.filename) == "png") {
            device <- function(..., width, height) {
                grDevices::png(..., width = 10, height = 10,
                               res = 300, units = "in")
            }
        } else if(tools::file_ext(input$survivalPlot.filename) == "pdf") {
            device <- function(..., width, height) {
                grDevices::pdf(..., width = 10, height = 10)
            }
        } else if(tools::file_ext(input$survivalPlot.filename) == "svg") {
            device <- function(..., width, height) {
                grDevices::svg(..., width = 10, height = 10)
            }
        }

        ggsave(file, plot = print(survival.plot(),newpage = F), device = device)
    })
