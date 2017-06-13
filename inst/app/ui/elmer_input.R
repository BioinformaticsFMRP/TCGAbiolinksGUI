tabItem(tabName = "elmerinput",
        fluidRow(
            column(8,  bsAlert("elmerinputmessage"),
                   bsCollapse(id = "collapseelmerinput", open = "Enhancer Linking by Methylation/Expression Relationship (ELMER) ",
                              bsCollapsePanel("Enhancer Linking by Methylation/Expression Relationship (ELMER) ",  includeHTML("elmer.html"), style = "default"),
                              bsCollapsePanel("Sample mapping matrix", dataTableOutput('elmerMaeSampleMapping'), style = "default"),
                              bsCollapsePanel("Sample metadata", dataTableOutput('elmerMaeSampleMetada'), style = "default"))
            ),
            column(4,
                   box(title = "Input - Create MAE object", width = NULL,
                       status = "danger",
                       solidHeader = FALSE, collapsible = FALSE,collapsed = FALSE,
                       box(title = "DNA methylation", width = NULL,
                           solidHeader = TRUE, collapsible = FALSE,collapsed = FALSE,
                           shinyFilesButton('elmermetfile', 'Select DNA methylation object', 'Please select DNA methylation object',
                                            multiple = FALSE),
                           bsTooltip("elmermetfile", "An R object (.rda) with DNA methylation data from HM450K platform (data frame or summarized experiment)",
                                     "left"),
                           selectizeInput('elmerInputMetPlatform',
                                          "Experiment group",
                                          choices=c("450K","EPIC"),
                                          multiple = FALSE),
                           numericInput("elmermetnacut", "DNA methylation: Cut-off NA samples (%)",
                                        min = 0, max = 1, value = 0.2, step = 0.1),
                           bsTooltip("elmermetnacut", "By default, for the DNA methylation data will remove probes with NA values in more than 20% samples and remove the anottation data.",
                                     "left")
                       ),  box(title = "Gene expression data", width = NULL,
                               solidHeader = TRUE, collapsible = FALSE,collapsed = FALSE,
                               shinyFilesButton('elmerexpfile', 'Select gene expression object', 'Please select gene expression object',
                                                multiple = FALSE),
                               bsTooltip("elmerexpfile", "An R object (.rda) with a gene expression object (data frame or summarized experiment)",
                                         "left"),
                               checkboxInput("elmerInputLinearizeExp", "Linearize gene expression values ?", value = FALSE, width = NULL),
                               bsTooltip("elmerInputLinearizeExp", "This will take the log2 of  gene expression to linearize the with DNA methylation","left")
                       ),  box(title = "Genome", width = NULL,
                               solidHeader = TRUE, collapsible = FALSE,collapsed = FALSE,
                               selectizeInput('elmerInputGenome',
                                              "Human reference genome",
                                              choices=c("GRCh38 (hg38)"="hg38","GRCh37 (hg19)"="hg19"),
                                              multiple = FALSE)
                       ),  box(title = "Annotation", width = NULL,
                               solidHeader = TRUE, collapsible = FALSE,collapsed = FALSE,
                               checkboxInput("elmerInputTCGA", "Is TCGA data ?", value = TRUE, width = NULL),
                               bsTooltip("elmerInputTCGA", "If checked we will automatically map each of the matrix (DNA methylation and gene expression) to the samples",
                                         "left"),
                               downloadButton('elmerInputSampleMapDownload', 'Download sample map example TSV'),
                               tags$br(),
                               tags$br(),
                               fileInput('elmerInputSampleMapFile', 'Select a tsv file with  sample mapping',
                                         accept = c(
                                             'text/csv',
                                             'text/comma-separated-values',
                                             'text/tab-separated-values',
                                             '.tsv',
                                             ".csv"
                                         )
                               ),
                               downloadButton('elmerInputpDataDownload', 'Download sample metadata example TSV'),
                               tags$br(),
                               tags$br(),
                               fileInput('elmerInputpDataFile', 'Select a tsv file with  sample metadata',
                                         accept = c(
                                             'text/csv',
                                             'text/comma-separated-values',
                                             'text/tab-separated-values',
                                             '.tsv',
                                             ".csv"
                                         )
                               )
                       ),
                       textInput("maesavefilename", "Save as:", value = "ELMER_input.rda", width = NULL, placeholder = NULL),
                       actionButton("elmercreatemae",
                                    "Create MAE object",
                                    style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                    icon = icon("floppy-o"))

                   )
            )
        )
)
