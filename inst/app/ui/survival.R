tabItem(tabName = "tcgasurvival",
        fluidRow(
            column(8,  bsAlert("survivalmessage"),
                   bsCollapse(id = "collapsesurvivalplot", open = "Survival plot",
                              bsCollapsePanel("Survival plot", uiOutput("survivalplot"), style = "default")
                   )),
            column(4,
                   box(title = "Survival plot",width = NULL,
                       status = "danger",
                       solidHeader = FALSE, collapsible = FALSE,
                       box(title = "Data",width = NULL,
                           solidHeader = TRUE, collapsible = TRUE,
                           bsTooltip("survivalplotfile", "A csv or rda file","left"),
                           shinyFilesButton('survivalplotfile', 'Select file', 'Please select a csv/rda file with a data frame with the following columns days_to_death, days_to_last_follow_up, vital_status',
                                            multiple = FALSE)),
                       box(title = "Parameters",width = NULL,
                           solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                           useShinyjs(),
                           checkboxInput("survivalbyGene", "Survival by  levels ?", value = FALSE, width = NULL),
                           numericInput("survivalPercent", "High/Low  cut-off (%)",
                                        min = 0, max = 1, value = 0.25, step = 0.1),
                           selectizeInput('survivalGene',
                                          'Gene/Probe',
                                          choices=NULL,
                                          multiple = FALSE),
                           selectizeInput('survivalplotgroup',
                                          'Group column',
                                          choices=NULL,
                                          multiple = FALSE),
                           textInput("survivalplotLegend", label = "Legend text", value = "Legend"),
                           textInput("survivalplotMain", label = "Title", value = "Kaplan-Meier Overall Survival Curves"),
                           bsTooltip("survivalplotLimit", "Set the limit of the x-axis, if 0 the automatic value will be considered",
                                     "left"),
                           useShinyjs(),
                           sliderInput("survivalplotLimit", "x-axis limit", min = 0, max = 10000, value = 0),
                           checkboxInput("survivalplotPvalue", "Add p-value?", value = TRUE, width = NULL),
                           checkboxInput("survivalplotRiskTable", "Add risk table?", value = TRUE, width = NULL),
                           checkboxInput("survivalplotConfInt", "Add confidence interval?", value = TRUE, width = NULL)
                       ),
                       box(title = "Size control",width = NULL,
                           solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                           sliderInput("survivalwidth", "Plot Width (%)", min = 0, max = 100, value = 100),
                           sliderInput("survivalheight", "Plot Height (px)", min = 0, max = 1200, value = 800)),
                       actionButton("survivalplotBt",
                                    "Plot survival plot",
                                    style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                    icon = icon("picture-o")))
            )
        )
)
