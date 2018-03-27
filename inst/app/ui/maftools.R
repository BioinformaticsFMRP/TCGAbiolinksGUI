tabItem(tabName = "maftools",
        fluidRow(
            column(8,  bsAlert("maftoolsmessage"),
                   bsCollapse(id = "collapsemaftools", open = "Plots",
                              bsCollapsePanel("Survival log rank test table", DT::dataTableOutput('maftoolsSurvivalTbl'), style = "default"),
                              bsCollapsePanel("Plots", uiOutput("maftoolsplot"), style = "default")
                   )),
            column(4,
                   box(title = "Maftools Plot",width = NULL,
                       status = "danger",
                       solidHeader = FALSE, collapsible = FALSE,
                       box(title = "Data",width = NULL,
                           solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                           shinyFilesButton('maffile2', 'Select MAF file', 'Please select a maf file',
                                            multiple = FALSE),
                           tags$br(),
                           tags$br(),
                           shinyFilesButton('clinicalmaf', 'Select Clinical file', 'Please select a clinical file',
                                            multiple = FALSE)
                       ),
                       box(title = "Plot selection",width = NULL,
                           solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                           selectizeInput('maftoolsPlot',
                                          "Plot selection",
                                          choices = c("plotmafSummary"="plotmafSummary",
                                                      "oncoplot"="oncoplot",
                                                      "Transition and Transversions"="titv",
                                                      "rainfallPlot"="rainfallPlot",
                                                      "survivalPlot"="survivalPlot",
                                                      "geneCloud"="geneCloud"),
                                          multiple = FALSE),
                           numericInput("maftoolsTop", "Top n mutated genes",
                                        min = 1, max = 300, value = 10, step = 1),
                           selectizeInput('maftoolsgene',
                                          "Specify gene",
                                          choices=NULL,
                                          multiple = FALSE),
                           selectizeInput('maftoolstsb',
                                          "Specify sample name",
                                          choices=NULL,
                                          multiple = FALSE),
                           numericInput("rainfallPlotps", "Point size",
                                        min = 0.1, max = 5, value = 0.6, step = 0.1),
                           checkboxInput("maftoolsrmNonMutated", "Remove non mutated", value = TRUE, width = NULL)
                       ),
                       actionButton("maftoolsPlotBt",
                                    "Plot",
                                    style = "background-color: #000080;
                                    color: #FFFFFF;
                                    margin-left: auto;
                                    margin-right: auto;
                                    width: 100%",
                                    icon = icon("picture-o")),
                       downloadButton('savemafpicture', 'Export figure', class = "butt2"),
                       textInput("mafPlot.filename", label = "Filename", value = "maftools.pdf"),
                       bsTooltip("mafPlot.filename", "Filename (pdf, png, svg)", "left"),
                       box(title = "Table selection",width = NULL,
                           solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                           numericInput("maftoolsTopMutated", "Evaluate survival for top n mutated genes",
                                        min = 1, max = 300, value = 10, step = 1),
                           actionButton("maftoolsTableBt",
                                        "Calculate",
                                        style = "background-color: #000080;
                                        color: #FFFFFF;
                                        margin-left: auto;
                                        margin-right: auto;
                                        width: 100%",
                                        icon = icon("picture-o"))
                       )

                   )
            )
        )
)
