tabItem(tabName = "gliomaclassifier",
        fluidRow(
            column(8, bsAlert("gliomaAlert"),
                   bsCollapse(id = "collapsegliomaClassifier", open = "Classification",
                              bsCollapsePanel("Classification", DT::dataTableOutput('gliomatbl'), style = "default")
                   )),
            column(4,
                   box(title = "Glioma classification",width = NULL,
                       status = "danger",
                       solidHeader = FALSE, collapsible = FALSE,
                       shinyFilesButton(id = 'classifyObj',
                                        label = 'Select DNA methylation (.rda) file',
                                        title = 'Please select a DNA methylation (.rda) file',
                                        multiple = FALSE),
                       tags$br(),
                       tags$br(),
                       actionButton("gliomaClassify",
                                    "Classify into gliomas molecular subtypes",
                                    style = "background-color: #000080;
                                    color: #FFFFFF;
                                    margin-left: auto;
                                    margin-right: auto;
                                    width: 100%",
                                    icon = icon("flask"))
                   )))
)
