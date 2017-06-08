tabItem(tabName = "tcgaSubtype",
        fluidRow(
            column(8, bsAlert("tcgaSubtypemessage"),
                   bsCollapse(id = "collapseTCGASubtype",
                              bsCollapsePanel("Subtype data: Summary",
                                              htmltools::div(style = "display:block;text-align: center;", plotlyOutput("subtypeview", width = 500, height = 500)),
                                              style = "default"),
                              bsCollapsePanel("Subtype data: results", dataTableOutput('tcgaSubtypetbl'), style = "default")
                   )),
            column(4,
                   box(title = "Subtype data search",width = NULL,
                       status = "danger",
                       solidHeader = FALSE, collapsible = FALSE,
                       selectizeInput('tcgasubtypeFilter',
                                      'Tumor filter',
                                      c("Adrenocortical carcinoma (TCGA-ACC)"="acc",
                                        "Breast invasive carcinoma (TCGA-BRCA)"="brca",
                                        "Colon adenocarcinoma (TCGA-COAD)"="coad",
                                        "Glioblastoma multiforme (TCGA-GBM)"="gbm",
                                        "Head and Neck squamous cell carcinoma (TCGA-HNSC)"="hnsc",
                                        "Kidney Chromophobe  (TCGA-KICH)"="kich",
                                        "Kidney renal papillary cell carcinoma (TCGA-KIRP)"="kirp",
                                        "Kidney renal clear cell carcinoma (TCGA-KIRC)"="kirc",
                                        "Brain Lower Grade Glioma (TCGA-LGG)"="lgg",
                                        "Lung adenocarcinoma  (TCGA-LUAD)"="luad",
                                        "Lung squamous cell carcinoma (TCGA-LUSC)"="lusc",
                                        "Prostate adenocarcinoma (TCGA-PRAD)"="prad",
                                        "Pheochromocytoma and Paraganglioma (TCGA-PCPG)"="pcpg",
                                        "pancan"="pancan",
                                        "Rectum adenocarcinoma (TCGA-READ)"="read",
                                        "Skin Cutaneous Melanoma (TCGA-SKCM)"="skcm",
                                        "Stomach adenocarcinoma (TCGA-STAD)"="stad",
                                        "Thyroid carcinoma (TCGA-THCA)"="thca",
                                        "Uterine Corpus Endometrial Carcinoma (TCGA-UCEC)"="ucec"),
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
                                    icon = icon("search")),
                       selectizeInput('subtypePlotCol',
                                      'Column to plot',
                                      NULL,
                                      multiple = FALSE)
                       )
            ))
)
