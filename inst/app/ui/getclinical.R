tabItem(tabName = "tcgaClinical",
        fluidRow(
            column(8, bsAlert("tcgaClinicalmessage"),
                   bsCollapse(id = "collapseTCGAClinical", open = "Description of the data",
                              bsCollapsePanel("Description of the data",  includeHTML("clinical_help.html"), style = "default"),
                              bsCollapsePanel("Clinical data", dataTableOutput('tcgaClinicaltbl'), style = "default")
                   )),
            column(4,
                   box(title = "Clinical data search",width = NULL,
                       status = "danger",
                       solidHeader = FALSE, collapsible = FALSE,
                       selectizeInput('tcgatumorClinicalFilter',
                                      'Tumor filter',
                                      NULL,
                                      multiple = FALSE),
                       checkboxInput("clinicalIndexed", "Indexed data ?", value = TRUE, width = NULL),
                       bsTooltip("clinicalIndexed", "The indexed data is a subset of the patient data updated with the last follow up. It does not have drug, radiation and other information, for that please use the XML data (not indexed)",
                                 "left"),
                       selectizeInput('tcgaClinicalTypeFilter',
                                      'Clinical file type filter',
                                      c("biospecimen", "clinical"),
                                      multiple = FALSE),
                       useShinyjs(),
                       selectizeInput('tcgaClinicalFilter',
                                      'Clinical information filter',
                                      c("drug", "admin", "follow_up", "radiation", "patient", "stage_event", "new_tumor_event"),
                                      multiple = FALSE),
                       bsTooltip("saveClinicalRda", "A Rda file stores a single R object (can only be openned in R)","left"),
                       checkboxInput("saveClinicalRda", "Save result as R object (Rdata) ?", value = FALSE, width = NULL),
                       checkboxInput("saveClinicalCsv", "Save result as csv ?", value = FALSE, width = NULL),
                       actionButton("tcgaClinicalBt",
                                    "Search",
                                    style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                    icon = icon("search")))
            )
        )
)
