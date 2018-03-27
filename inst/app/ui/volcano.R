tabItem(tabName = "volcano",
        fluidRow(
            column(8,  bsAlert("volcanomessage"),
                   bsCollapse(id = "collapseVolcano", open = "Volcano plot",
                              #bsCollapsePanel("Probes info", dataTableOutput('probesSE'), style = "default"),
                              bsCollapsePanel("Volcano plot",
                                              uiOutput("volcanoPlot"), style = "default")),
                   fluidRow(
                       column(width = 10,  offset = 1,
                              valueBoxOutput("volcanoBoxDown", width = 4),
                              valueBoxOutput("volcanoBoxInsig", width = 3),
                              valueBoxOutput("volcanoBoxUp", width = 4))
                   )),
            column(4,
                   box(title = "Volcano Plot",width = NULL,
                       status = "danger",
                       solidHeader = FALSE, collapsible = FALSE,
                       box(title = "Data",width = NULL,
                           solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                           bsTooltip("volcanofile", "A file with the DMR_ or DEA_ prefix. This is the differential expression analysis and differential methylation analysis outputs.", "left"),
                           shinyFilesButton('volcanofile', 'Select results', 'Please select the csv file with the results',
                                            multiple = FALSE)),
                       box(title = "Volcano options",width = NULL,
                           solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                           radioButtons("volcanoInputRb", "Select input type:",
                                        c("DNA methylation"="met",
                                          "Expression"="exp")),
                           numericInput("volcanoxcutMet", "DNA methylation threshold",
                                        min = 0, max = 1, value = 0, step = 0.05),
                           numericInput("volcanoxcutExp", "Log FC threshold",
                                        min = 0, max = 100, value = 0, step = 0.5),
                           numericInput("volcanoycut", "P-value adj cut-off",
                                        min = 0, max = 1, value = 0.05, step = 0.001)

                       ),
                       box(title = "Highlighting options",width = NULL,
                           solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                           checkboxInput("volcanoNames", "Show names?", value = FALSE, width = NULL),
                           checkboxInput("volcanoNamesFill", "Boxed names?", value = TRUE, width = NULL),
                           selectizeInput('volcanoHighlight',
                                          "Genes/Probes to Highlight",
                                          choices = NULL,  multiple = TRUE),
                           bsTooltip("volcanoShowHighlitgh", "Highlight genes/probes selected before, or significant probes/genes or both cases", "left"),
                           selectizeInput('volcanoShowHighlitgh',
                                          "Points to highlight",
                                          choices = c("highlighted","significant","both"),  multiple = FALSE)
                       ),
                       box(title = "Color control",width = NULL,
                           solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                           colourpicker::colourInput("volcanoColHighlight", "Highlight color", value = "orange"),
                           colourpicker::colourInput("colinsignificant", "Insignificant color", value = "black"),
                           colourpicker::colourInput("colHypomethylated", "Hypomethylated color", value = "darkgreen"),
                           colourpicker::colourInput("colHypermethylated", "Hypermethylate color", value = "red"),
                           colourpicker::colourInput("colDownregulated", "Downregulated color", value = "darkgreen"),
                           colourpicker::colourInput("colUpregulated", "Upregulated genes color", value = "red")
                       ),
                       box(title = "Size control",width = NULL,
                           solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                           sliderInput("volcanowidth", "Plot Width (%)", min = 0, max = 100, value = 100),
                           sliderInput("volcanoheight", "Plot Height (px)", min = 0, max = 1200, value = 400)),
                       bsTooltip("volcanoSave", "Along with the plot, it will save a new csv with the cut-offs set", "left"),
                       checkboxInput("volcanoSave", "Save file with results?", value = FALSE, width = NULL),
                       actionButton("volcanoPlotBt",
                                    "Volcano plot",
                                    style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                    icon = icon("picture-o")
                       ),
                       downloadButton('savevolvanopicture', 'Export figure', class = "butt2"),
                       textInput("volcanoPlot.filename", label = "Filename", value = "volcano.pdf"),
                       bsTooltip("volcanoPlot.filename", "Filename (pdf, png, svg)", "left")
                   )
            )
        )
)
