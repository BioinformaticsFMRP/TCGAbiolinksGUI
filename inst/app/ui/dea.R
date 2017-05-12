tabItem(tabName = "dea",
        fluidRow(
            column(8,  bsAlert("deamessage"),
                   bsCollapse(id = "collapsedea", open = "Genes info",
                              bsCollapsePanel("Genes info", dataTableOutput('deaSE'), style = "default")
                   )),
            column(4,
                   box(title = "DEA analysis",width = NULL,
                       status = "danger",
                       solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
                       box(title = "Data",width = NULL,
                           solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                           shinyFilesButton('deafile', 'Select SummarizedExperiment', 'Please select SummarizedExperiment object',
                                            multiple = FALSE)),
                       box(title = "Pre-analysis options",width = NULL,
                           solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                           checkboxInput("deanormalization", "Normalization of genes?", value = FALSE, width = NULL),
                           useShinyjs(),
                           selectizeInput('deanormalizationmet',
                                          "Normalization of genes method",
                                          choices = c("gcContent","geneLength"),
                                          multiple = FALSE),
                           checkboxInput("deafilter", "Quantile filter of genes?", value = FALSE, width = NULL),
                           selectizeInput('deafilteringmet',
                                          "DEA test method",
                                          choices = c("quantile","varFilter","filter1","filter2"),
                                          multiple = FALSE),
                           numericInput("deafilteringcut", "Threshold selected as mean for filtering",
                                        min = 0, max = 1, value = 0.25, step = 0.1)
                       ),
                       box(title = "Analysis parameter",width = NULL,
                           solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                           numericInput("deathrsld", "Log FC threshold",
                                        min = 0, max = 100, value = 0, step = 0.5),
                           numericInput("deapvalue", "P-value adj cut-off",
                                        min = 0, max = 1, value = 0.05, step = 0.001),
                           selectizeInput('deagroupCol',
                                          "Group column",
                                          choices = NULL,  multiple = FALSE),
                           selectizeInput('deagroup1',
                                          "Group 1",
                                          choices = NULL,  multiple = FALSE),
                           selectizeInput('deagroup2',
                                          "Group 2",
                                          choices = NULL,  multiple = FALSE),
                           selectizeInput('deamethod',
                                          "DEA test method",
                                          choices = c("glmLRT","exactTest"),
                                          multiple = FALSE),
                           actionButton("deaAnalysis",
                                        "dea analysis",
                                        style = "background-color: #000080;
                                        color: #FFFFFF;
                                        margin-left: auto;
                                        margin-right: auto;
                                        width: 100%",
                                        icon = icon("flask")))
                   )
            )
        )
)
