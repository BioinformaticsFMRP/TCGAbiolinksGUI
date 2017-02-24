tabItem(tabName = "pathview",
        fluidRow(
            column(8,  bsAlert("deamessage"),
                   bsCollapse(id = "collapsedea", open = "DEA plots",
                              bsCollapsePanel("Genes info", dataTableOutput('deaSE'), style = "default")
                   )),
            column(4,
                   box(title = "Pathway graphs",width = NULL,
                       status = "danger",
                       solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                       shinyFilesButton('pathewayexpfile', 'DEA result', 'Please select expression result object',
                                        multiple = FALSE),
                       selectizeInput('pathway.id',
                                      "pathway ID",
                                      choices = pathways.id,
                                      multiple = FALSE),
                       checkboxInput("kegg.native.checkbt", "Native KEGG?", value = TRUE, width = NULL),
                       sliderInput("pathwaygraphwidth", "Plot Width", min = 0, max = 1600, value = 1200),
                       sliderInput("pathwaygraphheight", "Plot Height", min = 0, max = 1200, value = 800),
                       actionButton("pathwaygraphBt",
                                    "Create pathway file",
                                    style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                    icon = icon("file-pdf-o")),
                       bsTooltip("pathwaygraphBt", "In the figure, colors represents the values of logFC. Red colors are genes up regulated in group 2 while green are the ones downregulated (for a file DEA_results_col_group1_group2...) ",
                                 "left")
                   )
            )
        )
)
