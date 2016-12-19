tabItem(tabName = "tcgaSubtype",
        fluidRow(
            column(8, bsAlert("tcgaSubtypemessage"),
                   bsCollapse(id = "collapseTCGASubtype", open = "Subtype data",
                              bsCollapsePanel("Subtype data", dataTableOutput('tcgaSubtypetbl'), style = "default")
                   )),
            column(4,
                   box(title = "Subtype data search",width = NULL,
                       status = "danger",
                       solidHeader = FALSE, collapsible = FALSE,
                       selectizeInput('tcgasubtypeFilter',
                                      'Tumor filter',
                                      c("Adrenocortical carcinoma"="acc",
                                        "Breast invasive carcinoma"="brca",
                                        "Colon adenocarcinoma"="coad",
                                        "Glioblastoma multiforme"="gbm",
                                        "Head and Neck squamous cell carcinoma"="hnsc",
                                        "Kidney Chromophobe"="kich",
                                        "Kidney renal papillary cell carcinoma"="kirp",
                                        "Kidney renal clear cell carcinoma"="kirc",
                                        "Brain Lower Grade Glioma"="lgg",
                                        "Lung adenocarcinoma"="luad",
                                        "Lung squamous cell carcinoma"="lusc",
                                        "Prostate adenocarcinoma"="prad",
                                        "pancan"="pancan",
                                        "Rectum adenocarcinoma"="read",
                                        "Skin Cutaneous Melanoma"="skcm",
                                        "Stomach adenocarcinoma"="stad",
                                        "Thyroid carcinoma"="thca",
                                        "Uterine Corpus Endometrial Carcinoma"="ucec"),
                                      multiple = FALSE),
                       bsTooltip("saveSubtypeRda", "A Rda file stores a single R object (can only be openned in R)","left"),
                       checkboxInput("saveSubtypeRda", "Save result as rda?", value = FALSE, width = NULL),
                       checkboxInput("saveSubtypeCsv", "Save result as csv?", value = FALSE, width = NULL),
                       actionButton("tcgaSubtypeBt",
                                    "Search",
                                    style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                    icon = icon("search")))
            ))
)
