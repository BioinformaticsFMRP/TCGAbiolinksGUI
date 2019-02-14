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
                       textInput("idatfilename", label = "Filename", value = "normalizedIdatSE_450K_hg38.rda"),
                       bsTooltip("idatfilename", "Extensions: csv (flat table) or rda (summarized Experiment)",
                                 "left"),
                       selectizeInput('IDATplatform',
                                      'DNA methylation platform',
                                      c( "450K" = "450K",
                                         "27K"  = "27K",
                                         "EPIC" = "EPIC"
                                      ),
                                      multiple = FALSE),
                       selectizeInput('IDATgenome',
                                      'Genome of reference to annotate probes',
                                      c( "hg38"  = "hg38",
                                         "hg19" = "hg19"
                                      ),
                                      multiple = FALSE),
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
