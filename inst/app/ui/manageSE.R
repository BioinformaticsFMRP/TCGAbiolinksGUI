tabItem(tabName = "seedit",
        fluidRow(
            column(8, bsAlert("seeditAlert"),
                   bsCollapse(id = "collapseseeEdit", open = "Description of the data",
                              bsCollapsePanel("Sample matrix", DT::dataTableOutput('seeditsampletb'), style = "default"),
                              bsCollapsePanel("Feature matrix", DT::dataTableOutput('seeditfeaturetb'), style = "default"),
                              bsCollapsePanel("Assay matrix", DT::dataTableOutput('seeditassaytb'), style = "default")
                   )),
            column(4,
                   box(title = "Editing summarized experiment",width = NULL,
                       status = "danger",
                       solidHeader = FALSE, collapsible = FALSE,
                       box(title = "Data",width = NULL,
                           solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                           shinyFilesButton('seeditfile', 'Select Summarized Experiment file', 'Please select summarized experiment object (.rda)', multiple = FALSE)
                       ),
                       box(title = "Update sample matrix",width = NULL,
                           solidHeader = TRUE, collapsible = TRUE,
                           downloadButton('downloadData', 'Download sample data matrix'),
                           tags$br(),
                           tags$br(),
                           fileInput('seuploadfile', 'Select a csv to update sample data matrix in the SE',
                                     accept = c(
                                         'text/csv',
                                         'text/comma-separated-values',
                                         'text/tab-separated-values',
                                         'text/plain',
                                         '.csv',
                                         '.rda'
                                     )
                           ),
                           textInput("seeditfilename", "Save as:", value = "edited_se.rda", width = "100%", placeholder = "Name of the new SummarizedExperiment file (.rda)"),
                           actionButton("seeditSEbt",
                                        "Save new Summarized Experiment (.rda) with the uploaded file",
                                        style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                        icon = icon("floppy-o")))
                   )))
)
