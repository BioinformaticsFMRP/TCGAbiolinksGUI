library(shiny)
library(shinyFiles)
library(TCGAbiolinks)
library(shinyBS)
library(shinyjs)
library(SummarizedExperiment)
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
        menuItem("biOMICs search" , tabName = "ontology", icon = icon("search")),
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
                    column(10, bsAlert("alert"),dataTableOutput('ontSearchtbl')),
                    column(2,
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
                                            icon = icon("search"))),
                           box(title = "Download",width = NULL,
                               status = "warning",
                               solidHeader = FALSE, collapsible = FALSE,
                               selectizeInput('ontftypeFilter',
                                              'Encode/Roadmap file filter',
                                              choices=NULL,
                                              multiple = TRUE),

                               shinyDirButton('folder', 'Download folder', 'Please select a folder',
                                              class='shinyDirectories btn-default', buttonType='warning'),

                               actionButton("ontSearchDownloadBt",
                                            "Download",
                                            style = "background-color: #F39C12;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 45%",
                                            icon = icon("download")),
                               bsTooltip("ontSearchDownloadBt", "Only downloads the files that were selected",
                                         "left"),
                               textOutput('ontSearchLink')
                           ),  box(title = "Report",width = NULL,
                                   status = "warning",
                                   solidHeader = FALSE, collapsible = FALSE,
                                   shinyDirButton('reportfolder', 'Report folder', 'Please select a folder where the report will be created',
                                                  class='shinyDirectories btn-default', buttonType='warning'),
                                   actionButton("ontReport",
                                                "Create report",
                                                style = "background-color: #F39C12;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 53%",
                                                icon = icon("book"))),
                           bsAlert("ontdownloaddirmessage")
                    )

                )
        ),
        tabItem(tabName = "tcgaSearch",
                fluidRow(
                    column(10, bsAlert("tcgasearchmessage"), dataTableOutput('tcgaSearchtbl')),
                    column(2,
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
                                              multiple = TRUE, selected = NULL),
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
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               selectizeInput('tcgasamplestypeFilter',
                                              'Sample type filter',
                                              table.code,
                                              multiple = TRUE),
                               conditionalPanel(
                                   condition = "input.tcgaExpFilter.indexOf('IlluminaHiSeq_RNASeqV2') > -1",
                                   selectizeInput('tcgaFrnaseqv2typeFilter',
                                                  'RNASeqV2 File type filter',
                                                  c("junction_quantification",
                                                    "rsem.genes.results",
                                                    "rsem.isoforms.results",
                                                    "rsem.genes.normalized_results",
                                                    "rsem.isoforms.normalized_results",
                                                    "bt.exon_quantification"),
                                                  multiple = TRUE)
                               ),  conditionalPanel(
                                   condition = "input.tcgaExpFilter.indexOf('IlluminaHiSeq_RNASeq') > -1",
                                   selectizeInput('tcgaFrnaseqtypeFilter',
                                                  'RNASeq File type filter',
                                                  c("exon.quantification",
                                                    "spljxn.quantification",
                                                    "gene.quantification"),
                                                  multiple = TRUE)
                               ),
                               conditionalPanel(
                                   condition = "input.tcgaExpFilter.indexOf('genome_wide_snp_6') > -1 && !input.tcgaExpFilter ",
                                   selectizeInput('tcgaFgwstypeFilter',
                                                  'genome_wide_snp_6 File type filter',
                                                  c("hg18.seg","hg19.seg","nocnv_hg18.seg","nocnv_hg19.seg"),
                                                  multiple = TRUE)
                               ),

                               shinyDirButton('tcgafolder', 'Folder select', 'Please select a folder',
                                              class='shinyDirectories btn-default', buttonType='warning'),
                               actionButton("tcgaDownloadBt",
                                            "Download",
                                            style = "background-color: #F39C12;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 50%",
                                            icon = icon("download"))
                           ),

                           box(title = "Prepare",width = NULL,
                               status = "warning",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               radioButtons("prepareRb", "Data type:",
                                            c("SummarizedExperiment" = TRUE,
                                              "Dataframe" = FALSE)),
                               conditionalPanel(
                                   condition = "input.prepareRb == 'TRUE'",
                                   checkboxInput("addSubTypeTCGA", "Add information about subtypes", value = FALSE, width = NULL)
                               ),
                               textInput("tcgafilename", "File name", value = "data.rda", width = NULL, placeholder = NULL),
                               shinyDirButton('tcgapreparefolder', 'Folder to save', 'Please select a folder',
                                              class='shinyDirectories btn-default', buttonType='warning'),
                               actionButton("tcgaPrepareBt",
                                            "Prepare data",
                                            style = "background-color: #F39C12;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 50%",
                                            icon = icon("cogs"))),
                           bsAlert("tcgaddirmessage")
                    )
                )
        ),
        tabItem(tabName = "tcgaOncoPrint",
                fluidRow(
                    column(10,  bsAlert("oncomessage"),
                           bsCollapse(id = "collapseOnco", open = "MAF files to download",
                                      bsCollapsePanel("MAF files to download", dataTableOutput('maftbl'), style = "default"),
                                      bsCollapsePanel("Oncoprint", uiOutput("oncoPlot"), style = "default")
                           )),
                    column(2,
                           box(title = "Download",width = NULL,
                               status = "warning",
                               solidHeader = FALSE, collapsible = FALSE,
                               shinyDirButton('maffolder', 'Folder select', 'Please select a folder',
                                              class='shinyDirectories btn-default', buttonType='warning'),
                               actionButton("mafDownloadBt",
                                            "Download",
                                            style = "background-color: #F39C12;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 50%",
                                            icon = icon("download"))
                           ),
                           bsAlert("oncoddirmessage"),
                           box(title = "Oncoprint",width = NULL,
                               status = "warning",
                               solidHeader = FALSE, collapsible = FALSE,
                               #fileInput('maffile', 'Choose maf File',
                               #          accept=c(".maf")),
                               shinyFilesButton('maffile', 'Select maf file', 'Please select a maf file',
                                                multiple = FALSE, buttonType='warning'),
                               selectizeInput('oncoGenes',
                                              "genes",
                                              choices = NULL,  multiple = TRUE),
                               colourInput("colDEL", "DEL colour", value = "red"),
                               colourInput("colINS", "INS colour", value = "blue"),
                               colourInput("colSNP", "SNP colour", value = "green"),
                               colourInput("colDNP", "DNP colour", value = "purple"),
                               sliderInput("oncowidth", "Plot Width (%)", min = 0, max = 100, value = 100),
                               sliderInput("oncoheight", "Plot Height (px)", min = 0, max = 800, value = 400),
                               actionButton("oncoprintPlot",
                                            "Plot oncoprint",
                                            style = "background-color: #F39C12;
                                              color: #FFFFFF;
                                              margin-left: auto;
                                              margin-right: auto;
                                              width: 100%",
                                            icon = icon("picture-o"))
                           )
                    )

                )
        ),
        tabItem(tabName = "dmr",

                fluidRow(
                    column(10,  bsAlert("dmrmessage"),
                           bsCollapse(id = "collapseDmr", open = "DMR plots",
                                      bsCollapsePanel("Probes info", dataTableOutput('probesSE'), style = "default"),
                                      bsCollapsePanel("DMR plots", uiOutput("dmrPlot"), style = "default"))),
                    column(2,
                           box(title = "DNA methylation object",width = NULL,
                               status = "warning",
                               solidHeader = FALSE, collapsible = FALSE,
                               shinyFilesButton('dmrfile', 'Select SummarizedExperiment', 'Please select SummarizedExperiment object',
                                                multiple = FALSE, buttonType='warning')),
                           #fileInput('dmrfile', 'Select SummarizedExperiment object',
                           #         accept=c(".rda"))),
                           box(title = "DMR analysis",width = NULL,
                               status = "warning",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
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
                               actionButton("volcanoPlot",
                                            "Volcano plot",
                                            style = "background-color: #F39C12;
                                              color: #FFFFFF;
                                              margin-left: auto;
                                              margin-right: auto;
                                              width: 100%",
                                            icon = icon("eye"))
                           ),
                           box(title = "Mean DNA methylation",width = NULL,
                               status = "warning",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               selectizeInput('meanmetgroupCol',
                                              "Group column",
                                              choices = NULL,  multiple = FALSE),
                               selectizeInput('meanmetsubgroupCol',
                                              "Sub group column",
                                              choices = NULL,  multiple = FALSE),
                               sliderInput("meanmetwidth", "Plot Width (%)", min = 0, max = 100, value = 100),
                               sliderInput("meanmetheight", "Plot Height (px)", min = 0, max = 800, value = 400),
                               actionButton("meanmetPlot",
                                            "DNA mean methylation plot",
                                            style = "background-color: #F39C12;
                                              color: #FFFFFF;
                                              margin-left: auto;
                                              margin-right: auto;
                                              width: 100%",
                                            icon = icon("eye"))

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
