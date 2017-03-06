tabItem(tabName = "netinf",
        fluidRow(
          column(8,  bsAlert("networkmessage"),
              bsCollapse(id = "collapseTCGAnetwork", open = "Top-scoring edges",
                              bsCollapsePanel("Top-scoring edges", dataTableOutput('tcgaNetworktbl'), style = "default")
                   )),
            column(4,
                   box(title = "Network inference",width = NULL,
                       status = "danger",
                       solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
                       box(title = "Data",width = NULL,
                           solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                           shinyFilesButton('networkfile', 'Select file',
                                            'Please select the summarized experiment file',
                                            multiple = FALSE)
                           ),
                        box(title = "Inference algorithm",width = NULL,
                           solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,

                           radioButtons("networkInferenceMethodRb", NULL,
                                        c("ARACNE"="aracne",
                                          "C3NET"="c3net",
                                           "CLR"= "clr",
                                           "MRNET"="mrnet")
                          )

                        ),
                        actionButton("networkBt",
                                    "Network inference",
                                    style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                    icon = icon("picture-o")
                       )
                      )
            )
        )
)
