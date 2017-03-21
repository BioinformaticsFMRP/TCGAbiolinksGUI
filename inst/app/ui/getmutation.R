tabItem(tabName = "tcgaMutation",
        fluidRow(
            column(8, bsAlert("tcgaMutationmessage"),
                   bsCollapse(id = "collapseTCGAMutation", open = "Mutation data",
                              bsCollapsePanel("Mutation data",div(dataTableOutput("tcgaMutationtbl"), style = "font-size: 90%;"), style = "default")
                   )),
            column(4,
                   box(title = "Mutatation data search",width = NULL,
                       status = "danger",
                       solidHeader = FALSE, collapsible = FALSE,
                       selectizeInput('tcgaMafTumorFilter',
                                      'Tumor filter',
                                      NULL,
                                      multiple = FALSE),
                       selectizeInput('tcgaMafPipeline',
                                      'Variant calling pipelines',
                                      choices = c("muse"="muse",
                                                  "somaticsniper"="somaticsniper",
                                                  "mutect"="mutect",
                                                  "varscan2"="varscan2"),
                                      multiple = FALSE),
                       checkboxInput("saveMafcsv", "Save MAF as csv?", value = TRUE, width = NULL),
                       actionButton("tcgaMafSearchBt",
                                    "Download",
                                    style = "background-color: #000080;
                                    color: #FFFFFF;
                                    margin-left: auto;
                                    margin-right: auto;
                                    width: 100%",
                                    icon = icon("download")))
            ))
)
