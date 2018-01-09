tabItem(tabName = "dnametidat",
        fluidRow(
            column(8, bsAlert("idatAlert"),
                   bsCollapse(id = "collapseIdat", open = "Samples identified",
                              bsCollapsePanel("Samples identified", DT::dataTableOutput('idattbl'), style = "default")
                   )),
            column(4,
                   box(title = "IDAT normalization",width = NULL,
                       status = "danger",
                       solidHeader = FALSE, collapsible = FALSE,
                       shinyDirButton('IDATfolder', 'Select idat folder', 'Please select a folder', TRUE),
                       tags$hr(),
                       textInput("idatfilename", label = "Filename", value = "normalizedIdatSE.rda"),
                       bsTooltip("idatfilename", "Extensions: csv (flat table) or rda (summarized Experiment)",
                                 "left"),
                       actionButton("idatqc",
                                    "Create shinyMethylSet object for QC visualization",
                                    style = "background-color: #000080;
                                    color: #FFFFFF;
                                    margin-left: auto;
                                    margin-right: auto;
                                    width: 100%",
                                    icon = icon("flask")),
                       bsTooltip("idatqc", "Use with shinyMethyl to visualizing the data and qc plots (function: runShinyMethyl)  ",
                                 "left"),
                       actionButton("idatnormalize",
                                    "Normalize and save",
                                    style = "background-color: #000080;
                     color: #FFFFFF;
                     margin-left: auto;
                     margin-right: auto;
                     width: 100%",
                                    icon = icon("flask")),
                       tags$hr(),
                       downloadButton('idatdownloadDataBut', 'Download normalized results',
                                      style = "background-color: #000080;
                       color: #FFFFFF;
                       margin-left: auto;
                       margin-right: auto;
                       width: 100%")
                   )
            )
        )
)
