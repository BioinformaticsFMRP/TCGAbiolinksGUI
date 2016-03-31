library(shiny)
library(shinyFiles)
library(TCGAbiolinks)
library(shinyBS)
library(shinyjs)
library(SummarizedExperiment)
library(pathview)
data(paths.hsa)
pathways.id <- names(paths.hsa)
names(pathways.id) <- unname(paths.hsa)

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

inputTextarea <- function(inputId, value="", nrows, ncols) {
    tagList(
        singleton(tags$head(tags$script(src = "textarea.js"))),
        tags$textarea(id = inputId,
                      class = "inputtextarea",
                      rows = nrows,
                      cols = ncols,
                      as.character(value))
    )
}


header <- dashboardHeader(
    title = "TCGAbiolinks"
)

sidebar <-  dashboardSidebar(
    sidebarMenu(
        #menuItem("biOMICs search" , tabName = "ontology", icon = icon("search")),
        #menuItem("Report" , tabName = "report", icon = icon("book")),
        menuItem("TCGA search" , tabName = "tcgaSearch", icon = icon("search")),
        menuItem("OncoPrint" , tabName = "tcgaOncoPrint", icon = icon("picture-o")),
        menuItem("Profile plot" , tabName = "tcgaProfilePlot", icon = icon("picture-o")),
        menuItem("Survival plot" , tabName = "tcgasurvival", icon = icon("picture-o")),
        menuItem("DMR analysis" , tabName = "dmr", icon = icon("flask")),
        menuItem("DEA analysis" , tabName = "dea", icon = icon("flask")),
        menuItem("Starburst plot" , tabName = "starburst", icon = icon("picture-o")),
        menuItem("Enrichment analysis" , tabName = "ea", icon = icon("flask"))
        #menuItem("ELMER analysis" , tabName = "elmer", icon = icon("flask"))
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
                    column(10,  bsAlert("alert"),
                           bsCollapse(id = "collapseonto", open = "Search results",
                                      bsCollapsePanel("Search results", dataTableOutput('ontSearchtbl'), style = "default"),
                                      bsCollapsePanel("biomicssummary", uiOutput("biomicssummary"), style = "default")
                           )),
                    column(2,
                           box(title = "Advanced search",width = NULL,
                               status = "danger",
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
                                            style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                            icon = icon("search"))),
                           box(title = "Download",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = FALSE,
                               selectizeInput('ontftypeFilter',
                                              'Encode/Roadmap file filter',
                                              choices=NULL,
                                              multiple = TRUE),

                               shinyDirButton('folder', 'Download folder', 'Please select a folder',
                                              class='shinyDirectories btn-default'),

                               actionButton("ontSearchDownloadBt",
                                            "Download",
                                            style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 45%",
                                            icon = icon("download")),
                               bsTooltip("ontSearchDownloadBt", "Only downloads the files that were selected",
                                         "left"),
                               textOutput('ontSearchLink')
                           ),  box(title = "Report",width = NULL,
                                   status = "danger",
                                   solidHeader = FALSE, collapsible = FALSE,
                                   shinyDirButton('reportfolder', 'Report folder', 'Please select a folder where the report will be created',
                                                  class='shinyDirectories btn-default'),
                                   actionButton("ontReport",
                                                "Create report",
                                                style = "background-color: #000080;
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
                    column(10, bsAlert("tcgasearchmessage"),
                           bsCollapse(id = "collapseTCGA", open = "TCGA search results",
                                      bsCollapsePanel("TCGA search results", dataTableOutput('tcgaSearchtbl'), style = "default")
                           )),
                    column(2,
                           box(title = "Advanced search",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE,
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
                                            style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                            icon = icon("search"))),
                           box(title = "Download",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               radioButtons("tcgaDownloadTypeRb", "Filter by:",
                                            c("No filter"="none",
                                              "Samples type"="type",
                                              "Samples barcode"="barcode")),
                               selectizeInput('tcgasamplestypeFilter',
                                              'Sample type filter',
                                              table.code,
                                              multiple = TRUE),
                               bsTooltip("tcgaDownloadBarcode", "Barcodes separeted by (;), (,) or (new line)",
                                         "left"),
                               useShinyjs(),
                               inputTextarea('tcgaDownloadBarcode', '', 2, 35),
                               useShinyjs(),
                               selectizeInput('tcgaFrnaseqv2typeFilter',
                                              'RNASeqV2 File type filter',
                                              c("junction_quantification",
                                                "rsem.genes.results",
                                                "rsem.isoforms.results",
                                                "rsem.genes.normalized_results",
                                                "rsem.isoforms.normalized_results",
                                                "bt.exon_quantification"),
                                              multiple = TRUE),
                               useShinyjs(),
                               selectizeInput('tcgaFrnaseqtypeFilter',
                                              'RNASeq File type filter',
                                              c("exon.quantification",
                                                "spljxn.quantification",
                                                "gene.quantification"),
                                              multiple = TRUE),
                               useShinyjs(),
                               selectizeInput('tcgaFgwstypeFilter',
                                              'genome_wide_snp_6 File type filter',
                                              c("hg18.seg","hg19.seg","nocnv_hg18.seg","nocnv_hg19.seg"),
                                              multiple = TRUE),
                               actionButton("tcgaDownloadBt",
                                            "Download",
                                            style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                            icon = icon("download"))
                           ),

                           box(title = "Prepare",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               radioButtons("prepareRb", "Data type:",
                                            c("SummarizedExperiment" = TRUE,
                                              "Dataframe" = FALSE)),
                               useShinyjs(),
                               checkboxInput("addSubTypeTCGA", "Add information about subtypes", value = FALSE, width = NULL),
                               textInput("tcgafilename", "File name", value = "data.rda", width = NULL, placeholder = NULL),
                               actionButton("tcgaPrepareBt",
                                            "Prepare data",
                                            style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                            icon = icon("cogs"))),
                           box(title = "Subtype search",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               selectizeInput('tcgasubtypeFilter',
                                              'Tumor filter',
                                              c("brca"="brca",
                                                "coad"="coad",
                                                "gbm"="gbm",
                                                "hnsc"="hnsc",
                                                "kich"="kich",
                                                "kirp"="kirp",
                                                "kirc"="kirc",
                                                "lgg"="lgg",
                                                "luad"="luad",
                                                "lusc"="lusc",
                                                "prad"="prad",
                                                "pancan"="pancan",
                                                "read"="read",
                                                "skcm"="skcm",
                                                "stad"="stad",
                                                "thca"="thca", "ucec"="ucec"),
                                              multiple = FALSE),
                               checkboxInput("saveSubtype", "Save as rda?", value = FALSE, width = NULL),
                               actionButton("tcgaSubtypeBt",
                                            "TCGA Subtype Search",
                                            style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                            icon = icon("search"))),
                           box(title = "Clinical data search",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               radioButtons("clinicalSearchType", "Search by:",
                                            c("Tumor" = TRUE,
                                              "Samples" = FALSE)),
                               useShinyjs(),
                               selectizeInput('tcgatumorClinicalFilter',
                                              'Tumor filter',
                                              unique(TCGAquery()$Disease),
                                              multiple = FALSE),
                               bsTooltip("clinicalBarcode", "Barcodes separeted by (;), (,) or (new line)",
                                         "left"),
                               useShinyjs(),
                               inputTextarea('clinicalBarcode', '', 2, 35),
                               selectizeInput('tcgaClinicalFilter',
                                              'Clinical file type filter',
                                              c("biospecimen_aliquot",
                                                "biospecimen_analyte",
                                                "biospecimen_cqcf",
                                                "biospecimen_diagnostic_slides",
                                                "biospecimen_normal_control",
                                                "biospecimen_portion",
                                                "biospecimen_protocol",
                                                "biospecimen_sample",
                                                "biospecimen_shipment_portion",
                                                "biospecimen_slide",
                                                "biospecimen_tumor_sample",
                                                "clinical_cqcf",
                                                "clinical_drug",
                                                "clinical_follow_up_v1.5",
                                                "clinical_follow_up_v2.1",
                                                "clinical_follow_up_v4.0",
                                                "clinical_follow_up_v4.0_nte",
                                                "clinical_nte",
                                                "clinical_omf_v4.0",
                                                "clinical_patient",
                                                "clinical_radiation"),
                                              multiple = FALSE),
                               checkboxInput("saveClinical", "Save result as rda?", value = FALSE, width = NULL),
                               actionButton("tcgaClinicalBt",
                                            "TCGA Subtype Search",
                                            style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                            icon = icon("search"))),
                           box(title = "MAF data search",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               selectizeInput('tcgaMafTumorFilter',
                                              'Tumor filter',
                                              unique(TCGAquery()$Disease),
                                              multiple = FALSE),
                               actionButton("tcgaMafSearchBt",
                                            "Search",
                                            style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 49%",
                                            icon = icon("search")),
                               actionButton("mafDownloadBt",
                                            "Download",
                                            style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 49%",
                                            icon = icon("download"))
                           ),
                           box(title = "Directory to save files",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
                               bsTooltip("tcgafolder", "Select a folder where data will be saved/downloaded",
                                         "left"),
                               shinyDirButton('tcgafolder', 'Folder to save', 'Please select a folder where files will be saved',
                                              class='shinyDirectories btn-default')
                           ),
                           bsAlert("tcgaddirmessage")
                    ))
        ),
        tabItem(tabName = "tcgaOncoPrint",
                fluidRow(
                    column(10,  bsAlert("oncomessage"),
                           bsCollapse(id = "collapseOnco", open = "Oncoprint",
                                      bsCollapsePanel("Oncoprint", uiOutput("oncoPlot"), style = "default")
                           )),
                    column(2,
                           box(title = "Oncoprint data",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
                               shinyFilesButton('maffile', 'Select maf file', 'Please select a maf file',
                                                multiple = FALSE),
                               radioButtons("oncoInputRb", "Genes by:",
                                            c("Selection"="Selection",
                                              "Text"="text")),
                               bsTooltip("oncoGenesTextArea", "Genes separeted by (;), (,) or (new line)",
                                         "left"),
                               useShinyjs(),
                               inputTextarea('oncoGenesTextArea', '', 2, 35),
                               selectizeInput('oncoGenes',
                                              "genes",
                                              choices = NULL,  multiple = TRUE)
                           ),
                           box(title = "Oncoprint metadata",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
                               shinyFilesButton('mafAnnotation', 'Select annotation file', 'Please select a file with the annotation data frame ',
                                                multiple = FALSE),
                               useShinyjs(),
                               selectizeInput('mafAnnotationcols',
                                              "Annotation columns",
                                              choices = NULL,  multiple = TRUE),
                               selectizeInput('mafAnnotationpos',
                                              "Annotation position",
                                              choices = c("top","bottom"),selected = "top",  multiple = FALSE)
                           ),
                           box(title = "Colors control",width = NULL,  status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               colourInput("colDEL", "DEL color", value = "red"),
                               colourInput("colINS", "INS color", value = "blue"),
                               colourInput("colSNP", "SNP color", value = "green"),
                               colourInput("colDNP", "DNP color", value = "purple")),
                           box(title = "Size control",width = NULL,  status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               sliderInput("oncowidth", "Plot Width (%)", min = 0, max = 100, value = 100),
                               sliderInput("oncoheight", "Plot Height (px)", min = 0, max = 800, value = 400)
                           ),
                           box(title = "Oncoprint plot",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
                               bsTooltip("oncoRmCols", "If there is no alteration in that sample, whether remove it on the oncoprint",
                                         "left"),
                               checkboxInput("oncoRmCols", "Remove empty columns?", value = FALSE, width = NULL),
                               checkboxInput("oncoShowColsNames", "Show column names?", value = FALSE, width = NULL),
                               actionButton("oncoprintPlot",
                                            "Plot oncoprint",
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
        ),
        tabItem(tabName = "dmr",

                fluidRow(
                    column(10,  bsAlert("dmrmessage"),
                           bsCollapse(id = "collapseDmr", open = "DMR plots",
                                      bsCollapsePanel("Probes info", dataTableOutput('probesSE'), style = "default"),
                                      bsCollapsePanel("DMR plots", uiOutput("dmrPlot"), style = "default"))),
                    column(2,
                           box(title = "DNA methylation object",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = FALSE,
                               shinyFilesButton('dmrfile', 'Select SummarizedExperiment', 'Please select SummarizedExperiment object',
                                                multiple = FALSE)),
                           #fileInput('dmrfile', 'Select SummarizedExperiment object',
                           #         accept=c(".rda"))),
                           box(title = "DMR analysis",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               sliderInput("dmrcores", "Cores",step=1,
                                           min = 1, max = parallel::detectCores(), value = 1),
                               numericInput("dmrthrsld", "DNA methylation threshold",
                                            min = 0, max = 1, value = 0, step = 0.05),
                               numericInput("dmrpvalue", "P-value adj cut-off",
                                            min = 0, max = 1, value = 0.05, step = 0.001),

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
                                            style = "background-color: #000080;
                                              color: #FFFFFF;
                                              margin-left: auto;
                                              margin-right: auto;
                                              width: 100%",
                                            icon = icon("flask"))),
                           box(title = "Volcano plot",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               colourInput("colHypomethylated", "Hypomethylated color", value = "darkgreen"),
                               colourInput("colHypermethylated", "Hypermethylate color", value = "red"),
                               colourInput("colinsignificant", "Insignificant color", value = "black"),
                               checkboxInput("dmrNamesVolcano", "Add probe names?", value = FALSE, width = NULL),
                               checkboxInput("dmrNamesVolcanoFill", "Fill names?", value = TRUE, width = NULL),
                               actionButton("volcanoPlot",
                                            "Volcano plot",
                                            style = "background-color: #000080;
                                              color: #FFFFFF;
                                              margin-left: auto;
                                              margin-right: auto;
                                              width: 100%",
                                            icon = icon("eye"))
                           ),
                           box(title = "Mean DNA methylation",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               selectizeInput('meanmetgroupCol',
                                              "Group column",
                                              choices = NULL,  multiple = FALSE),
                               selectizeInput('meanmetsubgroupCol',
                                              "Sub group column",
                                              choices = NULL,  multiple = FALSE),
                               checkboxInput("meanmetplotjitter", "Plot jitters?", value = TRUE, width = NULL),
                               conditionalPanel(
                                   condition = "input.meanmetplotjitter == TRUE",
                                   selectizeInput('meanmetsort',
                                                  "Sort method",
                                                  choices = c("None"=NULL,
                                                              "Ascending by mean"="mean.asc",
                                                              "Descending by mean"="mean.desc",
                                                              "Ascending by median"="median.asc",
                                                              "Descending by median"="median.desc"),
                                                  multiple = FALSE)
                               ),
                               sliderInput("meanmetAxisAngle", "Decimal:",   min = 0, max = 360, value = 90, step= 45),
                               actionButton("meanmetPlot",
                                            "DNA mean methylation plot",
                                            style = "background-color: #000080;
                                              color: #FFFFFF;
                                              margin-left: auto;
                                              margin-right: auto;
                                              width: 100%",
                                            icon = icon("eye"))

                           ),
                           box(title = "Heatmap",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               selectizeInput('colmetadataheatmap',
                                              "Annotations to the columns",
                                              choices = NULL,  multiple = TRUE),
                               selectizeInput('rowmetadataheatmap',
                                              "Annotations for rows",
                                              choices = NULL,  multiple = TRUE),
                               checkboxInput("heatmap.clusterrows", "Cluster rows?", value = FALSE, width = NULL),
                               checkboxInput("heatmap.clustercol", "Cluster columns?", value = FALSE, width = NULL),
                               checkboxInput("heatmap.show.row.names", "Show row names?", value = FALSE, width = NULL),
                               checkboxInput("heatmap.show.col.names", "Show col names?", value = FALSE, width = NULL),
                               actionButton("heatmapPlot",
                                            "Heatmap plot",
                                            style = "background-color: #000080;
                                              color: #FFFFFF;
                                              margin-left: auto;
                                              margin-right: auto;
                                              width: 100%",
                                            icon = icon("eye"))
                           ),
                           box(title = "Plot controls",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               sliderInput("meanmetwidth", "Plot Width (%)", min = 0, max = 100, value = 100),
                               sliderInput("meanmetheight", "Plot Height (px)", min = 0, max = 800, value = 400))
                    )
                )
        ),
        tabItem(tabName = "ea",
                fluidRow(
                    column(10,
                           bsAlert("eamessage"),
                           bsCollapse(id = "collapseEA", open = "EA plots",
                                      bsCollapsePanel("EA plots", uiOutput("eaPlot"), style = "default"))),
                    column(2,
                           box(title = "EA analysis",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = FALSE,
                               radioButtons("tcgaEaInputRb", "Genes by:",
                                            c("Selection"="Selection",
                                              "Text"="text")),
                               bsTooltip("eaGenesTextArea", "Genes separeted by (;), (,) or (new line)",
                                         "left"),
                               useShinyjs(),
                               inputTextarea('eaGenesTextArea', '', 2, 35),
                               selectizeInput('eagenes',
                                              "Genes",
                                              choices = unique(rownames(TCGAbiolinks:::EAGenes)),
                                              multiple = TRUE)),
                           box(title = "Colors control",width = NULL,  status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               colourInput("colBP", "Biological Process", value = "orange"),
                               colourInput("colCC", "Cellular Component color", value = "cyan"),
                               colourInput("colMF", "Molecular function color", value = "green"),
                               colourInput("colPat", "Pathways color", value = "yellow")),
                           box(title = "Plot control",width = NULL,  status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               sliderInput("nBar", "Number of bar histogram to show",
                                           step=1, min = 1, max = 20, value = 10),
                               sliderInput("eawidth", "Plot Width (%)", min = 0, max = 100, value = 100),
                               sliderInput("eaheight", "Plot Height (px)", min = 0, max = 1000, value = 1000)

                           ),
                           actionButton("eaplot",
                                        "EA barplot",
                                        style = "background-color: #000080;
                                              color: #FFFFFF;
                                              margin-left: auto;
                                              margin-right: auto;
                                              width: 100%",
                                        icon = icon("eye"))
                    )
                )
        ),
        tabItem(tabName = "tcgaProfilePlot",
                fluidRow(
                    column(10,  bsAlert("profileplotmessage"),
                           bsCollapse(id = "collapseprofileplot", open = "Profile plot",
                                      bsCollapsePanel("Profile plot", uiOutput("profileplot"), style = "default")
                           )),
                    column(2,
                           box(title = "Profile plot",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = FALSE,
                               shinyFilesButton('profileplotfile', 'Select rda file', 'Please select a rda file with a data frame',
                                                multiple = FALSE),
                               useShinyjs(),
                               selectizeInput('profileplotgroup',
                                              'Column with the group information',
                                              choices=NULL,
                                              multiple = FALSE),
                               selectizeInput('profileplotsubtype',
                                              'Column with the subtype information',
                                              choices=NULL,
                                              multiple = FALSE),
                               checkboxInput("profileplotrmnagroup", " Remove the NA groups?", value = FALSE, width = NULL),
                               checkboxInput("profileplotrmnasub", " Remove the NA subtypes?", value = FALSE, width = NULL)
                           ),
                           box(title = "Plot controls",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               numericInput("margin1", "Move upper line of vertical bar",
                                            min = -10, max = 10, value = 0.0, step = 0.1),
                               numericInput("margin2", "Move right line of vertical bar",
                                            min = -10, max = 10, value = 0.0, step = 0.1),
                               numericInput("margin3", "Move bottom line of vertical bar",
                                            min = -10, max = 10, value = 0.0, step = 0.1),
                               numericInput("margin4", "Move left line of vertical bar",
                                            min = -10, max = 10, value = 0.0, step = 0.1),
                               sliderInput("profilewidth", "Plot Width (%)", min = 0, max = 100, value = 100),
                               sliderInput("profileheight", "Plot Height (px)", min = 0, max = 800, value = 800)),
                           actionButton("profileplotBt",
                                        "Plot profile plot",
                                        style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                        icon = icon("eye"))
                    )
                )
        ),
        tabItem(tabName = "tcgasurvival",
                fluidRow(
                    column(10,  bsAlert("survivalmessage"),
                           bsCollapse(id = "collapsesurvivalplot", open = "survival plot",
                                      bsCollapsePanel("survival plot", uiOutput("survivalplot"), style = "default")
                           )),
                    column(2,
                           box(title = "survival plot",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = FALSE,
                               shinyFilesButton('survivalplotfile', 'Select rda file', 'Please select a rda file with a data frame',
                                                multiple = FALSE),
                               selectizeInput('survivalplotgroup',
                                              'Column with the group information',
                                              choices=NULL,
                                              multiple = FALSE)

                           ),
                           box(title = "Plot controls",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               sliderInput("survivalwidth", "Plot Width (%)", min = 0, max = 100, value = 100),
                               sliderInput("survivalheight", "Plot Height (px)", min = 0, max = 800, value = 800)),
                           actionButton("survivalplotBt",
                                        "Plot survival plot",
                                        style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                        icon = icon("eye"))
                    )
                )
        ),
        tabItem(tabName = "dea",

                fluidRow(
                    column(10,  bsAlert("deamessage"),
                           bsCollapse(id = "collapsedea", open = "DEA plots",
                                      bsCollapsePanel("Genes info", dataTableOutput('deaSE'), style = "default"),
                                      bsCollapsePanel("DEA plots", uiOutput("deaPlot"), style = "default"))),
                    column(2,
                           box(title = "Gene expression object",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = FALSE,
                               shinyFilesButton('deafile', 'Select SummarizedExperiment', 'Please select SummarizedExperiment object',
                                                multiple = FALSE)),
                           box(title = "Normalization of genes",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               checkboxInput("deanormalization", "Normalization of genes?", value = FALSE, width = NULL),
                               useShinyjs(),
                               selectizeInput('deanormalizationmet',
                                              "Normalization of genes method",
                                              choices = c("gcContent","geneLength"),
                                              multiple = FALSE)
                           ),
                           box(title = "Quantile filter of genes",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               checkboxInput("deafilter", "Quantile filter of genes?", value = FALSE, width = NULL),

                               selectizeInput('deafilteringmet',
                                              "DEA test method",
                                              choices = c("quantile","varFilter","filter1","filter2"),
                                              multiple = FALSE),
                               numericInput("deafilteringcut", "Threshold selected as mean for filtering",
                                            min = 0, max = 1, value = 0.25, step = 0.1)
                           ),
                           box(title = "dea analysis",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               numericInput("deathrsld", "Log FC threshold",
                                            min = 0, max = 1, value = 0, step = 0.05),
                               numericInput("deapvalue", "P-value adj cut-off",
                                            min = 0, max = 1, value = 0.05, step = 0.001),
                               selectizeInput('deagroupCol',
                                              "Group column",
                                              choices = NULL,  multiple = FALSE),
                               selectizeInput('deagroup1',
                                              "Group 1",
                                              choices = NULL,  multiple = FALSE),
                               selectizeInput('deagroup2',
                                              "Group 2",
                                              choices = NULL,  multiple = FALSE),
                               selectizeInput('deamethod',
                                              "DEA test method",
                                              choices = c("glmLRT","exactTest"),
                                              multiple = FALSE),
                               actionButton("deaAnalysis",
                                            "dea analysis",
                                            style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                            icon = icon("flask"))),
                           box(title = "Volcano plot",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               colourInput("colUpregulated", "Upregulated genes color", value = "red"),
                               colourInput("colDownregulated", "Down regulated color", value = "darkgreen"),
                               colourInput("coldeainsignificant", "Insignificant color", value = "black"),
                               actionButton("volcanodeaPlot",
                                            "Volcano plot",
                                            style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                            icon = icon("eye"))
                           ),
                           box(title = "Pathway graphs",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               selectizeInput('pathway.id',
                                              "pathway ID",
                                              choices = pathways.id,
                                              multiple = FALSE),
                               checkboxInput("kegg.native.checkbt", "Native KEGG?", value = TRUE, width = NULL),
                               actionButton("pathwaygraphBt",
                                            "Create pathway file",
                                            style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                            icon = icon("file-pdf-o"))
                           ),
                           box(title = "Plot controls",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               sliderInput("deawidth", "Plot Width (%)", min = 0, max = 100, value = 100),
                               sliderInput("deaheight", "Plot Height (px)", min = 0, max = 800, value = 400))
                    )
                )
        ),
        tabItem(tabName = "starburst",
                fluidRow(
                    useShinyjs(),
                    column(10,  bsAlert("starburstmessage"),
                           bsCollapse(id = "collapsestarburst", open = "starburst plots",
                                      bsCollapsePanel("starburst result - probe gene pairs", dataTableOutput('starburstResult'), style = "default"),
                                      bsCollapsePanel("starburst plots", uiOutput("starburstPlot"), style = "default"))),
                    column(2,
                           box(title = "Gene expression object",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = FALSE,
                               shinyFilesButton('starburstmetfile', 'Select SummarizedExperiment', 'Please select SummarizedExperiment object',
                                                multiple = FALSE),
                               shinyFilesButton('starburstexpfile', 'Select expression result', 'Please select expression result object',
                                                multiple = FALSE)),
                           box(title = "starburst analysis",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               numericInput("starburstexpFC", "Log FC threshold",
                                            min = 0, max = 10, value = 0, step = 0.05),
                               numericInput("starburstexFDR", "Expression FDR cut-off",
                                            min = 0, max = 1, value = 0.05, step = 0.001),
                               numericInput("starburstmetdiff", "Mean DNA methylation difference threshold",
                                            min = 0, max = 1, value = 0, step = 0.05),
                               numericInput("starburstmetFDR", "Methylation FDR cut-off",
                                            min = 0, max = 1, value = 0.05, step = 0.001),
                               selectizeInput('starburstgroup1',
                                              "Group 1",
                                              choices = NULL,  multiple = FALSE),
                               selectizeInput('starburstgroup2',
                                              "Group 2",
                                              choices = NULL,  multiple = FALSE),
                               checkboxInput("starburstNames", "Add genes names?", value = FALSE, width = NULL),
                               checkboxInput("starburstNamesFill", "Fill names?", value = TRUE, width = NULL)
                           ),
                           box(title = "Colors control",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               colourInput("sbcolInsignigicant", "Insignificant", value = "black"),
                               colourInput("sbcolUpHypo", "Upregulated & Hypomethylated", value =  "#E69F00"),
                               colourInput("sbcolDownHypo", "Downregulated & Hypomethylated", value = "#56B4E9"),
                               colourInput("sbcolHypo", "Hypomethylated", value = "#009E73"),
                               colourInput("sbcolHyper", "Hypermethylated", value = "red"),
                               colourInput("sbcolUp", "Upregulated", value =  "#0072B2"),
                               colourInput("sbcolDown", "Downregulated", value = "#D55E00"),
                               colourInput("sbcolUpHyper", "Upregulated & Hypermethylated", value = "#CC79A7"),
                               colourInput("sbcolDownHyper", "Downregulated & Hypermethylated", value = "purple")),
                           box(title = "Plot controls",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               sliderInput("starburstwidth", "Plot Width (%)", min = 0, max = 100, value = 100),
                               sliderInput("starburstheight", "Plot Height (px)", min = 0, max = 800, value = 400)),
                           actionButton("starburstPlot",
                                        "starburst plot",
                                        style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                        icon = icon("eye"))
                    )
                )
        ),
        tabItem(tabName = "elmer",
                fluidRow(
                    column(10,  bsAlert("elmermessage"),
                           bsCollapse(id = "collapsestelmer", open = "starburst plots",
                                      bsCollapsePanel("ELMER result", dataTableOutput('elmerResult'), style = "default"),
                                      bsCollapsePanel("Elmer plots", uiOutput("elmerPlot"), style = "default"))),
                    column(2,
                           box(title = "Create mee object", width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = FALSE,
                               shinyFilesButton('elmermetfile', 'Select DNA methylation object', 'Please select DNA methylation object',
                                                multiple = FALSE),
                               shinyFilesButton('elmerexpfile', 'Select expression object', 'Please select gene expression object',
                                                multiple = FALSE),
                               numericInput("elmermetnacut", "cut-off NA samples (%)",
                                            min = 0, max = 1, value = 0.2, step = 0.1),
                               actionButton("elmerpreparemee",
                                            "Create mee object",
                                            style = "background-color: #000080;
                                        color: #FFFFFF;
                                        margin-left: auto;
                                        margin-right: auto;
                                        width: 100%",
                                            icon = icon("floppy-o")))
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
