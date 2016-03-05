library(shiny)
library(shinyFiles)

table.code <- c('01','02','03','04','05','06','07','08','09','10',
                '11','12','13','14','20','40','50','60','61')
names(table.code) <- c("Primary solid Tumor","Recurrent Solid Tumor",
                       "Primary Blood Derived Cancer - Peripheral Blood",
                       "Recurrent Blood Derived Cancer - Bone Marrow",
                       "Additional - New Primary",
                       "Metastatic","Additional Metastatic",
                       "Human Tumor Original Cells",
                       "Primary Blood Derived Cancer - Bone Marrow",
                       "Blood Derived Normal","Solid Tissue Normal",
                       "Buccal Cell Normal","EBV Immortalized Normal",
                       "Bone Marrow Normal","Control Analyte",
                       "Recurrent Blood Derived Cancer - Peripheral Blood",
                       "Cell Lines","Primary Xenograft Tissue",
                       "Cell Line Derived Xenograft Tissue")

header <- dashboardHeader(
    title = "biOMICs"
)

sidebar <-  dashboardSidebar(
    sidebarMenu(
        menuItem("Ontology search" , tabName = "ontology", icon = icon("search")),
        #menuItem("Report" , tabName = "report", icon = icon("book")),
        menuItem("TCGA search" , tabName = "tcgaSearch", icon = icon("search"))
        #menuItem("ROADMAP search" , tabName = "roadmap", icon = icon("search")),
        #menuItem("ENCOCDE search" , tabName = "encode", icon = icon("search"))
    )
)

body <-  dashboardBody(
    tagList(
        singleton(tags$head(tags$script(
            paste0('//cdnjs.cloudflare.com/ajax/libs/datatables/',
                   '1.10.11/js/jquery.dataTables.min.js'),
            type = 'text/javascript'))),
        singleton(tags$head(tags$link(href = paste0('//cdn.datatables.net/',
                                                    'tabletools/2.2.4/css/',
                                                    'dataTables.tableTools.css'),
                                      rel = 'stylesheet',type = 'text/css'))),
        singleton(tags$head(tags$link(href = paste0('//cdn.datatables.net',
                                                    '/1.10.11/css/jquery.',
                                                    'dataTables.css'),
                                      rel = 'stylesheet',type = 'text/css'))),
        singleton(tags$head(tags$script(src = paste0('//cdn.datatables.net/',
                                                     'tabletools/2.2.4/js/dat',
                                                     'aTables.tableTools.min.js'),
                                        type = 'text/javascript'))),
        singleton(tags$head(tags$link(rel = "stylesheet", type = "text/css",
                                      href = "biOMICS.css")))
    ),
    tabItems(
        #tabItem(tabName = "report",
        #        fluidRow(includeHTML("main.html"))
        #),
        tabItem(tabName = "ontology",

                fluidRow(
                    column(9, dataTableOutput('ontSearchtbl')),
                    column(3,
                           box(title = "Advanced search",width = NULL,
                               status = "warning",
                               solidHeader = FALSE, collapsible = FALSE,
                               selectizeInput('ontExpFilter',
                                              'Experiments filter',
                                              platforms$Standard,
                                              multiple = TRUE),
                               selectizeInput('ontSamplesFilter',
                                              'Term',
                                              c("",
                                                union(
                                                    union(
                                                        roadmap.db$Sample.Name,
                                                        encode.db$biosample
                                                    ),
                                                    TCGAbiolinks::TCGAquery()$Disease
                                                )),
                                              multiple = FALSE,options = list(create = TRUE),selected = NULL),

                               shinyDirButton('folder', 'Folder select', 'Please select a folder',
                                              class='shinyDirectories btn-default', buttonType='warning'),
                               verbatimTextOutput("directorypath"),
                               actionButton("ontSearchBt",
                                            "search",
                                            style = "background-color: #F39C12;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                            icon = icon("search")),
                               verbatimTextOutput("system"),
                               selectizeInput('ontftypeFilter',
                                              'Encode/Roadmap file filter',
                                              unique(encode.db.files$file_format),
                                              multiple = TRUE),
                               actionButton("ontSearchDownloadBt",
                                            "Download selected",
                                            style = "background-color: #F39C12;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                            icon = icon("download")),

                               textOutput('ontSearchLink')
                           )
                    )

                )
        ),
        tabItem(tabName = "tcgaSearch",
                fluidRow(
                    column(9, dataTableOutput('tcgaSearchtbl')),
                    column(3,
                           box(title = "Advanced search",width = NULL,
                               status = "warning",
                               solidHeader = FALSE, collapsible = FALSE,
                               selectizeInput('tcgaTumorFilter',
                                              'Tumor filter',
                                              unique(TCGAquery()$Disease),
                                              multiple = TRUE),
                               selectizeInput('tcgaExpFilter',
                                              'Platforms filter',
                                              unique(TCGAquery()$Platform),
                                              multiple = TRUE),
                               selectizeInput('tcgaLevelFilter',
                                              'Level filter',
                                              c(1:3),
                                              multiple = TRUE, selected = NULL),
                               actionButton("tcgaSearchBt",
                                            "TCGA Search",
                                            style = "background-color: #F39C12;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                            icon = icon("search"))),
                           box(title = "Download",width = NULL,
                               status = "warning",
                               solidHeader = FALSE, collapsible = FALSE,
                               selectizeInput('tcgasamplestypeFilter',
                                              'Sample type filter',
                                              table.code,
                                              multiple = TRUE),
                               selectizeInput('tcgaFrnaseqv2typeFilter',
                                              'RNASeqV2 File type filter',
                                              c("junction_quantification",
                                                 "rsem.genes.results",
                                                "rsem.isoforms.results",
                                                "rsem.genes.normalized_results",
                                                "rsem.isoforms.normalized_results",
                                                "bt.exon_quantification"),
                                              multiple = TRUE),
                               selectizeInput('tcgaFrnaseqtypeFilter',
                                              'RNASeq File type filter',
                                              c("exon.quantification",
                                                "spljxn.quantification",
                                                "gene.quantification"),
                                              multiple = TRUE),
                               selectizeInput('tcgaFgwstypeFilter',
                                              'genome_wide_snp_6 File type filter',
                                              c("hg18.seg","hg19.seg","nocnv_hg18.seg","nocnv_hg19.seg"),
                                              multiple = TRUE),
                               shinyDirButton('tcgafolder', 'Folder select', 'Please select a folder',
                                              class='shinyDirectories btn-default', buttonType='warning'),
                               verbatimTextOutput("tcgadirectorypath"),
                               actionButton("tcgaDownloadBt",
                                            "Download selected",
                                            style = "background-color: #F39C12;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                            icon = icon("download")),
                               textOutput('tcgaSearchLink')
                           )
                    )
                )
        )
    )
)

# @title  Client side
# @description Client side - Download data from roadmap project
# @keywords internal
biOMICsUI <- dashboardPage(
    skin = "blue",
    header,
    sidebar,
    body
)
