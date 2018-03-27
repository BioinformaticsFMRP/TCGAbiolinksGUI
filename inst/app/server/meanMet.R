##----------------------------------------------------------------------
#                             Mean Methylation
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

observe({
    updateSelectizeInput(session, 'meanmetsubgroupCol', choices = {
        data <- meandata()
        if(!is.null(data)){
            data <- colData(data)
            #data <- data[,which(apply(data,2,function(x) class(parse_guess((x)))) == "character")]
            as.character(colnames(data))
        }
    }, server = TRUE)
})
observe({
    updateSelectizeInput(session, 'meanmetgroupCol', choices = {
        data <- meandata()
        if(!is.null(data)){
            data <- colData(data)
            #data <- data[,which(apply(data,2,function(x) class(parse_guess((x)))) == "character")]
            as.character(colnames(data))
        }
    }, server = TRUE)
})
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Plot
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
meanmetPlot <- reactive({
    input$meanmetPlot
    closeAlert(session, "meanmetAlert")
    jitter <- isolate({input$meanmetplotjitter})
    sort <- NULL
    if(isolate({input$meanmetSortCB})) sort <- isolate({input$meanmetsort})
    angle <- isolate({input$meanmetAxisAngle})
    data <- meandata()

    y.limits <- NULL
    if(isolate({input$meanmetSetLimits})) y.limits <- c(isolate({input$meanmetylimlow}),isolate({input$meanmetylimup}))

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

    legend.position <- isolate({input$meanmetLegendPos})
    legend.title.position <- isolate({input$meanmetLegendTitlePos})
    legend.ncols <- isolate({input$meanmetLegendncols})
    add.axis.x.text <- isolate({input$meanmetXnames})
    withProgress(message = 'Creating plot',
                 detail = 'This may take a while...', value = 0, {
                     if(is.null(sort)){
                         TCGAvisualize_meanMethylation(data=data,
                                                       groupCol=group,
                                                       subgroupCol=subgroup,
                                                       filename = NULL,
                                                       y.limits = y.limits,
                                                       legend.position = legend.position,
                                                       legend.title.position = legend.title.position,
                                                       legend.ncols = legend.ncols,
                                                       add.axis.x.text = add.axis.x.text,
                                                       plot.jitter = jitter,
                                                       axis.text.x.angle = angle )
                     } else {
                         TCGAvisualize_meanMethylation(data=data,
                                                       groupCol=group,
                                                       subgroupCol=subgroup,
                                                       filename = NULL,
                                                       legend.position = legend.position,
                                                       legend.title.position = legend.title.position,
                                                       legend.ncols = legend.ncols,
                                                       add.axis.x.text = add.axis.x.text,
                                                       y.limits = y.limits,
                                                       plot.jitter = jitter,
                                                       axis.text.x.angle = angle,
                                                       sort=sort)
                     }
                 })

})
observeEvent(input$meanmetPlot, {
    output$mean.plotting <- renderPlot({
        meanmetPlot()
    })
})

observeEvent(input$meanmetPlot , {
    updateCollapse(session, "collapsemeanmet", open = "Mean DNA methylation plot")
    output$meanMetplot <- renderUI({
        plotOutput("mean.plotting", width = paste0(isolate({input$meanmetwidth}), "%"), height = isolate({input$meanmetheight}))
    })})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# File selection
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observe({
    shinyFileChoose(input,
                    'meanmetfile',
                    roots = get.volumes(input$workingDir),
                    session = session,
                    restrictions = system.file(package ='base'),filetypes = c('', 'rda'))
})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Data input
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
meandata <-  reactive({
    inFile <- input$meanmetfile
    if (is.null(inFile)) return(NULL)
    file  <- as.character(parseFilePaths(get.volumes(isolate({input$workingDir})), input$meanmetfile)$datapath)

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


output$savemeanmetpicture <- downloadHandler(
    filename = function(){input$meanmetPlot.filename},
    content = function(file) {
        if(tools::file_ext(input$meanmetPlot.filename) == "png") {
            device <- function(..., width, height) {
                grDevices::png(..., width = 10, height = 10,
                               res = 300, units = "in")
            }
        } else if(tools::file_ext(input$meanmetPlot.filename) == "pdf") {
            device <- function(..., width, height) {
                grDevices::pdf(..., width = 10, height = 10)
            }
        } else if(tools::file_ext(input$meanmetPlot.filename) == "svg") {
            device <- function(..., width, height) {
                grDevices::svg(..., width = 10, height = 10)
            }
        }

        ggsave(file, plot = meanmetPlot(), device = device)
    })
