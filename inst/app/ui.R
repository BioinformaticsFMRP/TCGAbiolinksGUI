library(shiny)
library(shinyFiles)
library(TCGAbiolinks)
library(shinyBS)
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
        menuItem("TCGA search" , tabName = "tcgaSearch", icon = icon("search")),
        menuSubItem("TCGA - OncoPrint" , tabName = "tcgaOncoPrint"),
        menuItem("DMR analysis" , tabName = "dmr", icon = icon("flask"))
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
                    column(9, bsAlert("alert"),dataTableOutput('ontSearchtbl')),
                    column(3,
                           box(title = "Advanced search",width = NULL,
                               status = "warning",
                               solidHeader = FALSE, collapsible = FALSE,
                               selectizeInput('ontExpFilter',
                                              'Experiments filter',
                                              NULL,
                                              multiple = TRUE),
                               selectizeInput('ontSamplesFilter',
                                              'Term',
                                              choices=NULL,
                                              multiple = FALSE,options = list(create = TRUE),selected = NULL),


                               actionButton("ontSearchBt",
                                            "search",
                                            style = "background-color: #F39C12;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                            icon = icon("search")),
                               verbatimTextOutput("system")),
                           box(title = "Download",width = NULL,
                               status = "warning",
                               solidHeader = FALSE, collapsible = FALSE,
                               shinyDirButton('folder', 'Folder select', 'Please select a folder',
                                              class='shinyDirectories btn-default', buttonType='warning'),
                               verbatimTextOutput("directorypath"),
                               selectizeInput('ontftypeFilter',
                                              'Encode/Roadmap file filter',
                                              choices=NULL,
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
                           ),  box(title = "Report",width = NULL,
                                   status = "warning",
                                   solidHeader = FALSE, collapsible = FALSE,
                                   shinyDirButton('reportfolder', 'Folder for save the report', 'Please select a folder',
                                                  class='shinyDirectories btn-default', buttonType='warning'),
                                   verbatimTextOutput("reportdirectorypath"),
                                   actionButton("ontReport",
                                                "Create report",
                                                style = "background-color: #F39C12;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 49%",
                                                icon = icon("book"))),
                           singleton(
                               tags$head(tags$script(src = "message-handler.js"))
                           )
                    )

                )
        ),
        tabItem(tabName = "tcgaSearch",
                fluidRow(
                    column(9, bsAlert("tcgasearchmessage"), dataTableOutput('tcgaSearchtbl')),
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
                           ),
                           box(title = "Prepare",width = NULL,
                               status = "warning",
                               solidHeader = FALSE, collapsible = FALSE,
                               radioButtons("prepareRb", "Data type:",
                                            c("SummarizedExperiment" = TRUE,
                                              "Dataframe" = FALSE)),
                               shinyDirButton('tcgapreparefolder', 'Folder to save', 'Please select a folder',
                                              class='shinyDirectories btn-default', buttonType='warning'),
                               verbatimTextOutput("tcgapreparedir"),
                               actionButton("tcgaPrepareBt",
                                            "Prepare data",
                                            style = "background-color: #F39C12;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                            icon = icon("cogs")),
                               textOutput('tcgaprepare'))
                    )
                )
        ),
        tabItem(tabName = "tcgaOncoPrint",
                fluidRow(
                    column(9, dataTableOutput('maftbl'),
                           plotOutput("oncoPlot", click = "plot_click",height=800)),
                    column(3,
                           box(title = "Download",width = NULL,
                               status = "warning",
                               solidHeader = FALSE, collapsible = FALSE,
                               shinyDirButton('maffolder', 'Folder select', 'Please select a folder',
                                              class='shinyDirectories btn-default', buttonType='warning'),
                               verbatimTextOutput("mafdirectorypath"),
                               actionButton("mafDownloadBt",
                                            "Download selected",
                                            style = "background-color: #F39C12;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                            icon = icon("download")),
                               textOutput('mafDownloadfTxt')
                           ),box(title = "Oncoprint",width = NULL,
                                 status = "warning",
                                 solidHeader = FALSE, collapsible = FALSE,
                                 fileInput('maffile', 'Choose maf File',
                                           accept=c(".maf")),
                                 actionButton("oncoprintPlot",
                                              "Plot oncoprint",
                                              style = "background-color: #F39C12;
                                              color: #FFFFFF;
                                              margin-left: auto;
                                              margin-right: auto;
                                              width: 100%",
                                              icon = icon("picture-o")),
                                 selectizeInput('oncoGenes',
                                                "genes",
                                                choices = NULL,  multiple = TRUE)
                           )
                    )

                )
        ),
        tabItem(tabName = "dmr",

                fluidRow(
                    column(9, plotOutput("dmrPlot",height=800)),
                    column(3,
                           box(title = "DNA methylation object",width = NULL,
                               status = "warning",
                               solidHeader = FALSE, collapsible = FALSE,
                               fileInput('dmrfile', 'Select SummarizedExperiment object',
                                         accept=c(".rda"))),
                           box(title = "DMR analysis",width = NULL,
                               status = "warning",
                               solidHeader = FALSE, collapsible = FALSE,
                               numericInput("dmrthrsld", "DNA methylation threshold",
                                            min = 0, max = 1, value = 0, step = 0.05),
                               numericInput("dmrpvalue", "P-value adj cut-off",
                                            min = 0, max = 1, value = 0.05, step = 0.001),
                               sliderInput("dmrcores", "Cores",step=1,
                                           min = 1, max = parallel::detectCores(), value = 1),
                               selectizeInput('dmrgroupCol',
                                              "Group column",
                                              choices = NULL,  multiple = FALSE),
                               selectizeInput('dmrgroup1',
                                              "Group 1",
                                              choices = NULL,  multiple = FALSE),
                               selectizeInput('dmrgroup2',
                                              "Group 2",
                                              choices = NULL,  multiple = FALSE),
                               actionButton("dmrAnalysis",
                                            "DMR analysis",
                                            style = "background-color: #F39C12;
                                              color: #FFFFFF;
                                              margin-left: auto;
                                              margin-right: auto;
                                              width: 100%",
                                            icon = icon("flask")),
                               textOutput('tcgadmr')
                           ),
                           box(title = "Mean DNA methylation",width = NULL,
                               status = "warning",
                               solidHeader = FALSE, collapsible = FALSE,
                               selectizeInput('meanmetgroupCol',
                                              "Group column",
                                              choices = NULL,  multiple = FALSE),
                               selectizeInput('meanmetsubgroupCol',
                                              "Sub group column",
                                              choices = NULL,  multiple = FALSE),
                               actionButton("meanmetPlot",
                                            "DNA mean methylation plot",
                                            style = "background-color: #F39C12;
                                              color: #FFFFFF;
                                              margin-left: auto;
                                              margin-right: auto;
                                              width: 100%",
                                            icon = icon("eye")),
                               textOutput('tcgameanmet')
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
