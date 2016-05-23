library(shiny)
library(shinyFiles)
library(TCGAbiolinks)
library(shinyBS)
library(shinyjs)
library(SummarizedExperiment)
library(pathview)
library(reshape2)
data(paths.hsa)
pathways.id <- names(paths.hsa)
names(pathways.id) <- unname(paths.hsa)

menu.icon <- "arrow-circle-right"


disease <- TCGAbiolinks:::disease.table$abbreviation
names(disease) <- TCGAbiolinks:::disease.table$name
disease <- sort(disease)
disease <- disease[disease %in% unique(TCGAquery()$Disease)]

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
    title = "TCGAbiolinks",
    titleWidth = 250
)
header$children[[2]]$children <-  tags$a(href='http://bioconductor.org/packages/TCGAbiolinksGUI/',
                                         tags$img(src='logo_no_text.png',height='50',width='220'))

sidebar <-  dashboardSidebar(
    width = 250,
    sidebarMenu(
        tags$hr(class="lineData"),
        menuItem("TCGA Data Information" , tabName = "tcgaSummary", icon = icon("info")),
        menuItem("Get TCGA data", icon = icon("database"),
                 menuSubItem("Molecular data" , tabName = "tcgaSearch", icon = icon("database")),
                 menuSubItem("Mutation data" , tabName = "tcgaMutation", icon = icon("database")),
                 menuSubItem("Clinical data" , tabName = "tcgaClinical", icon = icon("database")),
                 menuSubItem("Subtype data" , tabName = "tcgaSubtype", icon = icon("database"))
        ),
        tags$hr(class="lineAnalysis"),
        menuItem("Clinical analysis", icon = icon("flask"),
                 menuSubItem("Profile plot" , tabName = "tcgaProfilePlot", icon = icon("picture-o")),
                 menuSubItem("Survival plot" , tabName = "tcgasurvival", icon = icon("picture-o"))
        ),

        menuItem("Epigenetic analysis", icon = icon("flask"),
                 menuSubItem("Differential methylation analysis" , tabName = "dmr", icon = icon("flask")),
                 menuSubItem("Volcano plot" , tabName = "volcano", icon = icon("picture-o")),
                 menuSubItem("Heatmap plot" , tabName = "heatmap", icon = icon("picture-o")),
                 menuItem("Mean DNA methylation plot" , tabName = "meanmet", icon = icon("picture-o"))
        ),
        menuItem("Transcriptomic analysis", icon =  icon("flask"),
                 menuSubItem("Differential expression analysis" , tabName = "dea", icon = icon("flask")),
                 menuSubItem("Volcano plot" , tabName = "volcano", icon = icon("picture-o")),
                 menuSubItem("Heatmap plot" , tabName = "heatmap", icon = icon("picture-o")),
                 menuSubItem("Enrichment analysis" , tabName = "ea", icon = icon("flask"))
        ),
        menuItem("Genomic analysis", icon = icon("flask"),
                 menuSubItem("OncoPrint plot" , tabName = "tcgaOncoPrint", icon = icon("picture-o"))
        ),
        tags$hr(class="lineIntegrative"),
        #menuItem("Integrative analysis", icon = icon("flask"),
        menuItem("Starburst plot" , tabName = "starburst", icon = icon("picture-o")),
        menuItem("ELMER" , tabName = "elmer", icon = icon("flask")),
        #),
        #menuItem("Starburst plot" , tabName = "starburst", icon = icon("picture-o")),
        #menuItem("ELMER analysis", icon = icon("flask"),
        #         menuSubItem("Create ELMER object" , tabName = "elmer", icon = icon("database")),
        #         menuSubItem("Analysis" , tabName = "elmer", icon = icon("flask")),
        #         menuSubItem("Visualize Results" , tabName = "elmer", icon = icon("picture-o"))
        #),
        #tags$hr(class="lineConfig"),
        #menuItem("Configuration", tabName = "config", icon = icon("cogs")),
        tags$hr(class="lineDoc"),
        menuItem("Tutorial/Vignettes", icon = icon("book"),
                 menuSubItem("TCGAbiolinksGUI Manual" , href = "https://www.bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/tcgaBiolinks.html", icon = icon("external-link")),
                 menuSubItem("TCGAbiolinks Manual" , href = "https://www.bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/tcgaBiolinks.html", icon = icon("external-link")),
                 menuSubItem("ELMER Manual" , href = "https://www.bioconductor.org/packages/3.3/bioc/vignettes/ELMER/inst/doc/vignettes.pdf", icon = icon("external-link"))
    ),
    menuItem("References", icon = icon("file-text-o"),
             menuSubItem("TCGAbiolinks" , href = "https://doi.org/10.1093/nar/gkv1507", icon = icon("external-link")),
             menuSubItem("ELMER" , href = "https://doi.org/10.1186/s13059-015-0668-3", icon = icon("external-link"))
    ), div(id = "greetbox-outer",
    menuItem("welcome", tabName = "welcome",selected = T)))
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
                                      href = "TCGAbiolinksGUI.css")))
    ),

    tabItems(
        tabItem(tabName = "welcome",
                fluidRow(
                    column(1),
                    column(10,
                    includeHTML("index.html")
                    ))),
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
                                              disease,
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
                               bsTooltip("tcgaDownloadBarcode", "Barcodes separeted by (;), (,) or (new line)", "left"),
                               useShinyjs(),
                               inputTextarea('tcgaDownloadBarcode', '', 2, 30),
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

                           box(title = "Directory to save & prepare files",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
                               bsTooltip("tcgafolder", "Select a folder where data will be read/saved/downloaded",
                                         "left"),
                               shinyDirButton('tcgafolder', 'Select directory', 'Please select a folder where files will be saved',
                                              class='shinyDirectories btn-default')
                           ),
                           bsAlert("tcgaddirmessage")
                    ))
        ),
        tabItem(tabName = "tcgaClinical",
                fluidRow(
                    column(10, bsAlert("tcgaClinicalmessage"),
                           bsCollapse(id = "collapseTCGAClinical", open = "Clinical data",
                                      bsCollapsePanel("Clinical data", dataTableOutput('tcgaClinicaltbl'), style = "default")
                           )),
                    column(2,
                           box(title = "Search",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = FALSE,
                               radioButtons("clinicalSearchType", "Search by:",
                                            c("Tumor" = TRUE,
                                              "Samples" = FALSE)),
                               useShinyjs(),
                               selectizeInput('tcgatumorClinicalFilter',
                                              'Tumor filter',
                                              disease,
                                              multiple = FALSE),
                               bsTooltip("clinicalBarcode", "Barcodes separeted by (;), (,) or (new line)",
                                         "left"),
                               useShinyjs(),
                               inputTextarea('clinicalBarcode', '', 2, 30),
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
                               bsTooltip("saveClinicalRda", "A Rda file stores a single R object (can only be openned in R)","left"),
                               checkboxInput("saveClinicalRda", "Save result as Rda?", value = FALSE, width = NULL),
                               checkboxInput("saveClinicalCsv", "Save result as csv?", value = FALSE, width = NULL),
                               actionButton("tcgaClinicalBt",
                                            "Search",
                                            style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                            icon = icon("search")))
                    ))),
        tabItem(tabName = "tcgaMutation",
                fluidRow(
                    column(10, bsAlert("tcgaMutationmessage"),
                           bsCollapse(id = "collapseTCGAMutation", open = "Mutation data",
                                      bsCollapsePanel("Mutation data", dataTableOutput('tcgaMutationtbl'), style = "default")
                           )),
                    column(2,
                           box(title = "Search: Mutation data",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = FALSE,
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
                                     width: 100%",
                                            icon = icon("search")),
                               actionButton("mafDownloadBt",
                                            "Download",
                                            style = "background-color: #000080;
                                     color: #FFFFFF;
                                     margin-left: auto;
                                     margin-right: auto;
                                     width: 100%",
                                            icon = icon("download")))
                    ))),
        tabItem(tabName = "tcgaSubtype",
                fluidRow(
                    column(10, bsAlert("tcgaSubtypemessage"),
                           bsCollapse(id = "collapseTCGASubtype", open = "Subtype data",
                                      bsCollapsePanel("Subtype data", dataTableOutput('tcgaSubtypetbl'), style = "default")
                           )),
                    column(2,
                           box(title = "Subtype search",width = NULL,
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
                    ))),
        tabItem(tabName = "tcgaSummary",
                fluidRow(
                    column(10, bsAlert("tcgaSummaryMessage"),
                           bsCollapse(id = "collapseTCGAsummary", open = "Summary",
                                      bsCollapsePanel("Summary", uiOutput("tcgaSummary"), style = "default")
                           )),
                    column(2,
                           box(title = "Summary plot",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = FALSE,
                               box(title = "Type of plot",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                                   radioButtons("summaryInputRb", "Plot type:",
                                                c("Type of samples"="type",
                                                  "Platforms vs samples"="platsample"))),
                               box(title = "Parameters",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                                   selectizeInput('tcgaSummaryTumorFilter',
                                                  'Tumor',
                                                  unique(TCGAquery()$Disease),
                                                  multiple = TRUE),
                                   selectizeInput('tcgaSummaryExpFilter',
                                                  'Platforms',
                                                  unique(TCGAquery()$Platform),
                                                  multiple = TRUE, selected = NULL),
                                   selectizeInput('tcgaSummarySamplestypeFilter',
                                                  'Sample type',
                                                  table.code,
                                                  multiple = TRUE),
                                   selectizeInput('tcgaSummaryLevelFilter',
                                                  'Level filter',
                                                  c(1:3),
                                                  multiple = FALSE, selected = 1),
                                   checkboxInput("summaryAddBarCount", "Add counts label to barchart?", value = TRUE, width = NULL),
                                   checkboxInput("summaryShowSampleType", "Highlight sample type in barchart?", value = FALSE, width = NULL)
                               ),
                               box(title = "Colors control",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   colourInput("summarySetsBarColor", "Sets Bar color", value = "#56B4E9"),
                                   colourInput("summaryQueryColorTP", " Primary solid Tumor color", value = colorRampPalette(rainbow(20),interpolate = "spline")(20)[1]),
                                   colourInput("summaryQueryColorTR", " Recurrent Solid Tumor color", value = colorRampPalette(rainbow(20),interpolate = "spline")(20)[2]),
                                   colourInput("summaryQueryColorTB", " Primary Blood Derived Cancer - Peripheral Blood color", value = colorRampPalette(rainbow(20),interpolate = "spline")(20)[3]),
                                   colourInput("summaryQueryColorTRBM", " Recurrent Blood Derived Cancer - Bone Marrow color", value = colorRampPalette(rainbow(20),interpolate = "spline")(20)[4]),
                                   colourInput("summaryQueryColorTAP", " Additional - New Primary color", value = colorRampPalette(rainbow(20),interpolate = "spline")(20)[5]),
                                   colourInput("summaryQueryColorTM", " Metastatic color", value = colorRampPalette(rainbow(20),interpolate = "spline")(20)[5]),
                                   colourInput("summaryQueryColorTAM", " Additional Metastatic color", value = colorRampPalette(rainbow(20),interpolate = "spline")(20)[6]),
                                   colourInput("summaryQueryColorTHOC", " Human Tumor Original Cells color", value = colorRampPalette(rainbow(20),interpolate = "spline")(20)[7]),
                                   colourInput("summaryQueryColorTBM", " Primary Blood Derived Cancer - Bone Marrow color", value = colorRampPalette(rainbow(20),interpolate = "spline")(20)[8]),
                                   colourInput("summaryQueryColorNB", " Blood Derived Normal color", value = colorRampPalette(rainbow(20),interpolate = "spline")(20)[9]),
                                   colourInput("summaryQueryColorNT", " Solid Tissue Normal color", value = colorRampPalette(rainbow(20),interpolate = "spline")(20)[10]),
                                   colourInput("summaryQueryColorNBC", " Buccal Cell Normal color", value = colorRampPalette(rainbow(20),interpolate = "spline")(20)[11]),
                                   colourInput("summaryQueryColorNEBV", " EBV Immortalized Normal color", value = colorRampPalette(rainbow(20),interpolate = "spline")(20)[12]),
                                   colourInput("summaryQueryColorNBM", " Bone Marrow Normal color", value = colorRampPalette(rainbow(20),interpolate = "spline")(20)[13]),
                                   colourInput("summaryQueryColorCELLC", " Control Analyte color", value = colorRampPalette(rainbow(20),interpolate = "spline")(20)[14]),
                                   colourInput("summaryQueryColorTRB", " Recurrent Blood Derived Cancer - Peripheral Blood color", value = colorRampPalette(rainbow(20),interpolate = "spline")(20)[15]),
                                   colourInput("summaryQueryColorCELL", " Cell Lines color", value = colorRampPalette(rainbow(20),interpolate = "spline")(20)[16]),
                                   colourInput("summaryQueryColorXP", " Primary Xenograft Tissue color", value = colorRampPalette(rainbow(20),interpolate = "spline")(20)[17]),
                                   colourInput("summaryQueryColorXCL", " Cell Line Derived Xenograft Tissue color", value = colorRampPalette(rainbow(20),interpolate = "spline")(20)[18])
                                   ),

                               box(title = "Size controls",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   #sliderInput("summaryncol", "Number columns", min = 1, max = 4, value = 3),
                                   sliderInput("summarywidth", "Plot Width (%)", min = 0, max = 100, value = 100),
                                   sliderInput("summaryheight", "Plot Height (px)", min = 0, max = 1200, value = 800)),
                               actionButton("tcgaSummaryBt",
                                            "Plot summary",
                                            style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                            icon = icon("picture-o")))
                    )
                )
        ),
        tabItem(tabName = "tcgaOncoPrint",
                fluidRow(
                    column(10,  bsAlert("oncomessage"),
                           bsCollapse(id = "collapseOnco", open = "Oncoprint",
                                      bsCollapsePanel("Oncoprint", uiOutput("oncoPlot"), style = "default")
                           )),
                    column(2,
                           box(title = "Oncoprint Plot",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = FALSE,
                               box(title = "Data",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                                   shinyFilesButton('maffile', 'Select MAF file', 'Please select a maf file',
                                                    multiple = FALSE),
                                   tags$br(),
                                   tags$br(),
                                   bsTooltip("mafAnnotation", "The object must have a column with the ID of patients. We look for one of the following columns: bcr_patient_barcode or patient","left"),
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
                               box(title = "Gene selection",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   radioButtons("oncoInputRb", "Genes by:",
                                                c("Selection"="Selection",
                                                  "Text"="text",
                                                  "File"="file")),
                                   bsTooltip("oncoGenesTextArea", "Genes separeted by (;), (,) or (new line)",
                                             "left"),
                                   useShinyjs(),
                                   inputTextarea('oncoGenesTextArea', '', 2, 30),
                                   selectizeInput('oncoGenes',
                                                  "genes",
                                                  choices = NULL,  multiple = TRUE),
                                   shinyFilesButton('oncoGenesFiles', 'Select file with genes', 'Please select a file with genes',  multiple = FALSE)
                               ),
                               box(title = "General parameters",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   selectizeInput('oncoHeatmapLegendSide',
                                                  "Heatmap legend side",
                                                  choices = c("right","bottom"),  multiple = FALSE),
                                   selectizeInput('oncoAnnotationLegendSide',
                                                  "Annotation legend side",
                                                  choices = c("right","bottom"),  multiple = FALSE),
                                   bsTooltip("oncoRmCols", "If there is no alteration in that sample, whether remove it on the oncoprint",
                                             "left"),
                                   checkboxInput("oncoRowSort", "Sort by frequency of alterations decreasingly?", value = FALSE, width = NULL),
                                   checkboxInput("oncoRmCols", "Remove empty columns?", value = FALSE, width = NULL),
                                   checkboxInput("oncoShowColsNames", "Show column names?", value = FALSE, width = NULL),
                                   checkboxInput("oncoShowRowBarplot", "Show barplot annotation on rows?", value = TRUE, width = NULL),
                                   numericInput("oncoWSpace", "Distance between columns",
                                                min = 0, max = 2, value = 0.5, step = 0.1),
                                   numericInput("oncoHSpace", "Distance between rows",
                                                min = 0, max = 2, value = 0.5, step = 0.1)
                               ),
                               box(title = "Text parameters",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   numericInput("oncoTextRowSize", "Rows text size",
                                                min = 0, max = 20, value = 10.0, step = 1),
                                   numericInput("oncoTextLabelSize", "Labels text Size",
                                                min = 0, max = 20, value = 10.0, step = 1)),
                               box(title = "Colors control",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   colourInput("colBG", "Background color", value = "#CCCCCC"),
                                   colourInput("colDEL", "DEL color", value = "red"),
                                   colourInput("colINS", "INS color", value = "blue"),
                                   colourInput("colSNP", "SNP color", value = "#008000"),
                                   colourInput("colDNP", "DNP color", value = "purple")),
                               box(title = "Size control",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   sliderInput("oncowidth", "Plot Width (%)", min = 0, max = 100, value = 100),
                                   sliderInput("oncoheight", "Plot Height (px)", min = 0, max = 2000, value = 800)
                               ),

                               actionButton("oncoprintPlot",
                                            "Plot oncoprint",
                                            style = "background-color: #000080;
                                              color: #FFFFFF;
                                              margin-left: auto;
                                              margin-right: auto;
                                              width: 100%",
                                            icon = icon("picture-o")
                               ))
                    )
                )
        ),
        tabItem(tabName = "volcano",
                fluidRow(
                    column(10,  bsAlert("volcanomessage"),
                           bsCollapse(id = "collapseVolcano", open = "Volcano plot",
                                      #bsCollapsePanel("Probes info", dataTableOutput('probesSE'), style = "default"),
                                      bsCollapsePanel("Volcano plot", uiOutput("volcanoPlot"), style = "default")),
                           fluidRow(
                               column(width = 10,  offset = 1,
                                      valueBoxOutput("volcanoBoxDown", width = 4),
                                      valueBoxOutput("volcanoBoxInsig", width = 3),
                                      valueBoxOutput("volcanoBoxUp", width = 4))
                           )),
                    column(2,
                           box(title = "Volcano Plot",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = FALSE,
                               box(title = "Data",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                                   shinyFilesButton('volcanofile', 'Select results', 'Please select the csv file with the results',
                                                    multiple = FALSE)),
                               box(title = "Volcano options",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   radioButtons("volcanoInputRb", "Select input type:",
                                                c("DNA methylation"="met",
                                                  "Expression"="exp")),
                                   numericInput("volcanoxcutMet", "DNA methylation threshold",
                                                min = 0, max = 1, value = 0, step = 0.05),
                                   numericInput("volcanoxcutExp", "Log FC threshold",
                                                min = 0, max = 100, value = 0, step = 0.5),
                                   numericInput("volcanoycut", "P-value adj cut-off",
                                                min = 0, max = 1, value = 0.05, step = 0.001)

                               ),
                               box(title = "Highligthing options",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   checkboxInput("volcanoNames", "Show names?", value = FALSE, width = NULL),
                                   checkboxInput("volcanoNamesFill", "Boxed names?", value = TRUE, width = NULL),
                                   selectizeInput('volcanoHighlight',
                                                  "Genes/Probes to Highlight",
                                                  choices = NULL,  multiple = TRUE),
                                   selectizeInput('volcanoShowHighlitgh',
                                                  "Points to highlight",
                                                  choices = c("highlighted","significant","both"),  multiple = FALSE)
                               ),
                               box(title = "Color control",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   colourInput("volcanoColHighlight", "Highlight color", value = "orange"),
                                   colourInput("colinsignificant", "Insignificant color", value = "black"),
                                   colourInput("colHypomethylated", "Hypomethylated color", value = "darkgreen"),
                                   colourInput("colHypermethylated", "Hypermethylate color", value = "red"),
                                   colourInput("colDownregulated", "Downregulated color", value = "darkgreen"),
                                   colourInput("colUpregulated", "Upregulated genes color", value = "red")
                               ),
                               box(title = "Size control",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   sliderInput("volcanowidth", "Plot Width (%)", min = 0, max = 100, value = 100),
                                   sliderInput("volcanoheight", "Plot Height (px)", min = 0, max = 1200, value = 400)),
                               bsTooltip("volcanoSave", "Along with the plot, it will save a new csv with the cut-offs set", "left"),
                               checkboxInput("volcanoSave", "Save file with results?", value = FALSE, width = NULL),
                               actionButton("volcanoPlotBt",
                                            "Volcano plot",
                                            style = "background-color: #000080;
                                 color: #FFFFFF;
                                 margin-left: auto;
                                 margin-right: auto;
                                 width: 100%",
                                            icon = icon("picture-o")
                               ))
                    )
                )
        ),
        tabItem(tabName = "meanmet",

                fluidRow(
                    column(10,  bsAlert("meanmetmessage"),
                           bsCollapse(id = "collapsemeanmet", open = "Mean DNA methylation plot",
                                      bsCollapsePanel("Mean DNA methylation plot", uiOutput("meanMetplot"), style = "default"))),
                    column(2,
                           box(title = "Mean DNA methylation",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = FALSE,
                               box(title = "Data",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   shinyFilesButton('meanmetfile', 'Select file', 'Please select SummarizedExperiment object',
                                                    multiple = FALSE)),
                               box(title = "Parameters control",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   selectizeInput('meanmetgroupCol',
                                                  "Group column",
                                                  choices = NULL,  multiple = FALSE),
                                   selectizeInput('meanmetsubgroupCol',
                                                  "Sub group column",
                                                  choices = NULL,  multiple = FALSE),
                                   checkboxInput("meanmetplotjitter", "Plot jitters?", value = TRUE, width = NULL),
                                   checkboxInput("meanmetSetLimits", "Select y limits?", value = FALSE, width = NULL),
                                   numericInput("meanmetylimup", "Y upper limit",
                                                step=0.1, min = 0, max = 1, value = 1),
                                   numericInput("meanmetylimlow", "Y lower limit",
                                                step=0.1, min = 0, max = 1, value = 0),
                                   checkboxInput("meanmetSortCB", "Sort plot?", value = FALSE, width = NULL),
                                   selectizeInput('meanmetsort',
                                                  "Sort method",
                                                  choices = c("None"=NULL,
                                                              "Ascending by mean"="mean.asc",
                                                              "Descending by mean"="mean.desc",
                                                              "Ascending by median"="median.asc",
                                                              "Descending by median"="median.desc"),
                                                  multiple = FALSE),
                                   sliderInput("meanmetAxisAngle", "x-axis label angle:",   min = 0, max = 360, value = 90, step= 45)),
                               box(title = "Size control",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   sliderInput("meanmetwidth", "Plot Width (%)", min = 0, max = 100, value = 100),
                                   sliderInput("meanmetheight", "Plot Height (px)", min = 0, max = 1200, value = 800)),
                               actionButton("meanmetPlot",
                                            "DNA mean methylation plot",
                                            style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                            icon = icon("picture-o")))
                    )
                )
        ),
        tabItem(tabName = "heatmap",

                fluidRow(
                    column(10,  bsAlert("heatmapmessage"),
                           bsCollapse(id = "collapseHeatmap", open = "Heatmap",
                                      bsCollapsePanel("Heatmap", uiOutput("heatmapPlot"), style = "default"))),
                    column(2,
                           box(title = "Heatmap",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = FALSE,
                               box(title = "Data",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                                   shinyFilesButton('heatmapfile', 'Select file',
                                                    'Please select SummarizedExperiment or data frame object',
                                                    multiple = FALSE),
                                   tags$br(),
                                   tags$br(),
                                   shinyFilesButton('heatmapresultsfile', 'Select results', 'Please select object',
                                                    multiple = FALSE)),
                               box(title = "Type of heatmap",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,

                                   radioButtons("heatmapTypeInputRb", NULL,
                                                c("DNA methylation"="met",
                                                  "Gene expression"="exp"),selected = "met")
                               ),
                               box(title = "Genes/Probes selection",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,

                                   useShinyjs(),
                                   radioButtons("heatmapProbesInputRb", "Select probes by:",
                                                c("Status"="Status",
                                                  "Text"="text")),
                                   inputTextarea('heatmapProbesTextArea', '', 2, 30),
                                   checkboxInput("heatmap.hypoprobesCb", "Hypermethylatd probes", value = TRUE, width = NULL),
                                   checkboxInput("heatmap.hyperprobesCb", "Hypomethylated probes", value = TRUE, width = NULL),
                                   radioButtons("heatmapGenesInputRb", "Select genes by:",
                                                c("Status"="Status",
                                                  "Text"="text")),
                                   inputTextarea('heatmapGenesTextArea', '', 2, 30),
                                   checkboxInput("heatmap.upGenesCb", "Up regulated genes", value = TRUE, width = NULL),
                                   checkboxInput("heatmap.downGenewsCb", "Down regulated genes", value = TRUE, width = NULL)),
                               box(title = "Annotations options",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,

                                   selectizeInput('colmetadataheatmap',
                                                  "Columns annotations",
                                                  choices = NULL,  multiple = TRUE),
                                   checkboxInput("heatmap.sortCb", "Sort by column?", value = FALSE, width = NULL),
                                   selectizeInput('heatmapSortCol',
                                                  "Sort by columns",
                                                  choices = NULL,  multiple = FALSE),
                                   selectizeInput('rowmetadataheatmap',
                                                  "Rows annotations",
                                                  choices = NULL,  multiple = TRUE)),
                               box(title = "Other options",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,

                                   selectizeInput('heatmapScale',
                                                  "Scale data",
                                                  choices = c("none","row","col"),selected = "none",  multiple = FALSE),
                                   checkboxInput("heatmap.clusterrows", "Cluster rows?", value = FALSE, width = NULL),
                                   checkboxInput("heatmap.clustercol", "Cluster columns?", value = FALSE, width = NULL),
                                   checkboxInput("heatmap.show.row.names", "Show row names?", value = FALSE, width = NULL),
                                   checkboxInput("heatmap.show.col.names", "Show col names?", value = FALSE, width = NULL)),
                               box(title = "Size control",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   sliderInput("heatmapwidth", "Plot Width (%)", min = 0, max = 100, value = 100),
                                   sliderInput("heatmapheight", "Plot Height (px)", min = 0, max = 800, value = 1000)
                               ),
                               actionButton("heatmapPlotBt",
                                            "Heatmap plot",
                                            style = "background-color: #000080;
                                                 color: #FFFFFF;
                                                 margin-left: auto;
                                                 margin-right: auto;
                                                 width: 100%",
                                            icon = icon("picture-o")))
                    )
                )
        ),
        tabItem(tabName = "dmr",

                fluidRow(
                    column(10,  bsAlert("dmrmessage"),
                           bsCollapse(id = "collapseDmr", open = "DMR plots",
                                      bsCollapsePanel("Probes info", dataTableOutput('probesSE'), style = "default"))),
                    column(2,
                           box(title = "DMR analysis",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = FALSE,
                               box(title = "Data",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE,
                                   bsTooltip("dmrfile", "A summarized Experiment", "left"),
                                   shinyFilesButton('dmrfile', 'Select data (.rda)', 'Please select SummarizedExperiment object',
                                                    multiple = FALSE)),
                               #fileInput('dmrfile', 'Select SummarizedExperiment object',
                               #         accept=c(".rda"))),
                               box(title = "Parameters control",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   sliderInput("dmrcores", "Cores",step=1,
                                               min = 1, max = parallel::detectCores(), value = 1),
                                   numericInput("dmrthrsld", "DNA methylation threshold",
                                                min = 0, max = 1, value = 0, step = 0.05),
                                   numericInput("dmrpvalue", "P-value adj cut-off",
                                                min = 0, max = 1, value = 0.05, step = 0.001),

                                   selectizeInput('dmrgroupCol',
                                                  "Group column",
                                                  choices = NULL,  multiple = FALSE),
                                   selectizeInput('dmrgroups',
                                                  "Groups",
                                                  choices = NULL,  multiple = TRUE)),
                               actionButton("dmrAnalysis",
                                            "DMR analysis",
                                            style = "background-color: #000080;
                                              color: #FFFFFF;
                                              margin-left: auto;
                                              margin-right: auto;
                                              width: 100%",
                                            icon = icon("flask")))
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
                               box(title = "Gene selection",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE,collapsed = FALSE,
                                   radioButtons("tcgaEaInputRb", "Genes by:",
                                                c("Selection"="Selection",
                                                  "Text"="text",
                                                  "File"="file")),
                                   bsTooltip("eaGenesTextArea", "Genes separeted by (;), (,) or (new line)",
                                             "left"),
                                   useShinyjs(),
                                   inputTextarea('eaGenesTextArea', '', 2, 30),
                                   selectizeInput('eagenes',
                                                  "Genes",
                                                  choices = unique(rownames(TCGAbiolinks:::EAGenes)),
                                                  multiple = TRUE),
                                   bsTooltip("eaGenesFiles", "Expecting a column Gene_symbol or mRNA in the file (csv,txt,rda). If it is the DEA_results files, it will ignore the insignificant genes",
                                             "left"),
                                   shinyFilesButton('eaGenesFiles', 'Select file with genes', 'Please select a file with genes',  multiple = FALSE)),
                               box(title = "Parameters control",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   numericInput("eaSizeText", "Size of the text",
                                                step=0.1, min = 0.1, max = 3, value = 1),
                                   numericInput("eaxlim", "X upper limit",
                                                step=1, min = 0, max = 50, value = 0),
                                   sliderInput("nBar", "Number of bar histogram to show",
                                               step=1, min = 1, max = 20, value = 10)),
                               box(title = "Plot selection",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   checkboxInput("eaPlotBPCB", " Plot biological process?", value = TRUE, width = NULL),
                                   checkboxInput("eaPlotCCCB", " Plot Cellular Component?", value = TRUE, width = NULL),
                                   checkboxInput("eaPlotMFCB", " Plot Molecular function?", value = TRUE, width = NULL),
                                   checkboxInput("eaPlotPatCB", " Plot Pathways?", value = TRUE, width = NULL)),
                               box(title = "Colors control",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   colourInput("colBP", "Biological Process", value = "orange"),
                                   colourInput("colCC", "Cellular Component color", value = "cyan"),
                                   colourInput("colMF", "Molecular function color", value = "green"),
                                   colourInput("colPat", "Pathways color", value = "yellow")),
                               box(title = "Size control",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   sliderInput("eawidth", "Plot Width (%)", min = 0, max = 100, value = 100),
                                   sliderInput("eaheight", "Plot Height (px)", min = 0, max = 1000, value = 800)),
                               actionButton("eaplot",
                                            "EA barplot",
                                            style = "background-color: #000080;
                                              color: #FFFFFF;
                                              margin-left: auto;
                                              margin-right: auto;
                                              width: 100%",
                                            icon = icon("picture-o"))
                           )
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
                           box(title = "Profle Plot ",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = FALSE,
                               box(title = "Data ",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE,
                                   shinyFilesButton('profileplotfile', 'Select file', 'Please select a rda file with a data frame',
                                                    multiple = FALSE)),
                               box(title = "General parameters ",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   useShinyjs(),
                                   selectizeInput('profileplotgroup',
                                                  'Group column',
                                                  choices=NULL,
                                                  multiple = FALSE),
                                   selectizeInput('profileplotsubtype',
                                                  'Subgruoup column',
                                                  choices=NULL,
                                                  multiple = FALSE),
                                   checkboxInput("profileplotrmnagroup", " Remove the NA groups?", value = FALSE, width = NULL),
                                   checkboxInput("profileplotrmnasub", " Remove the NA subgroups?", value = FALSE, width = NULL)),
                               box(title = "Text parameters",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   numericInput("profileplot.axis.title.size", "Size title left barplot",
                                                min = 0, max = 10, value = 1.0, step = 0.1),
                                   numericInput("profileplot.axis.textsize", "Text size cluster names",
                                                min = 0, max = 10, value = 1.0, step = 0.1),
                                   numericInput("profileplot.legend.size", "Legend text size",
                                                min = 0, max = 10, value = 1.0, step = 0.1),
                                   numericInput("profileplot.legend.title.size", "Title legend size",
                                                min = 0, max = 10, value = 1.0, step = 0.1),
                                   numericInput("profileplot.geom.label.size", "Size % cluster distribution",
                                                min = 0, max = 10, value = 4.0, step = 0.1)
                               ),
                               box(title = "Colors control",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   colourInput("profileplotColorTextLeftBar", "Sets % text color", value = "black")),
                               box(title = "Left bar control",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   numericInput("margin1", "Move upper line of vertical bar",
                                                min = -10, max = 10, value = 0.0, step = 0.1),
                                   numericInput("margin2", "Move right line of vertical bar",
                                                min = -10, max = 10, value = 0.0, step = 0.1),
                                   numericInput("margin3", "Move bottom line of vertical bar",
                                                min = -10, max = 10, value = 0.0, step = 0.1),
                                   numericInput("margin4", "Move left line of vertical bar",
                                                min = -10, max = 10, value = 0.0, step = 0.1)),
                               box(title = "Size control",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   sliderInput("profilewidth", "Plot Width (%)", min = 0, max = 100, value = 100),
                                   sliderInput("profileheight", "Plot Height (px)", min = 0, max = 1200, value = 800)),
                               actionButton("profileplotBt",
                                            "Plot profile plot",
                                            style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                            icon = icon("picture-o")))
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
                           box(title = "Survival plot",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = FALSE,
                               box(title = "Data",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE,
                                   bsTooltip("survivalplotfile", "A csv or rda file","left"),
                                   shinyFilesButton('survivalplotfile', 'Select file', 'Please select a csv/rda file with a data frame',
                                                    multiple = FALSE)),
                               box(title = "Parameters",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   selectizeInput('survivalplotgroup',
                                                  'Group column',
                                                  choices=NULL,
                                                  multiple = FALSE),
                                   textInput("survivalplotLegend", label = "Legend text", value = "Legend"),
                                   textInput("survivalplotMain", label = "Title", value = "Kaplan-Meier Overall Survival Curves"),
                                   bsTooltip("survivalplotLimit", "Set the limit of the x-axis, if 0 the automatic value will be considered",
                                             "left"),
                                   useShinyjs(),
                                   sliderInput("survivalplotLimit", "x-axis limit", min = 0, max = 10000, value = 0),
                                   checkboxInput("survivalplotPvalue", "Add p-value?", value = TRUE, width = NULL)

                               ),
                               box(title = "Size control",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   sliderInput("survivalwidth", "Plot Width (%)", min = 0, max = 100, value = 100),
                                   sliderInput("survivalheight", "Plot Height (px)", min = 0, max = 1200, value = 800)),
                               actionButton("survivalplotBt",
                                            "Plot survival plot",
                                            style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                            icon = icon("picture-o")))
                    )
                )
        ),
        tabItem(tabName = "dea",
                fluidRow(
                    column(10,  bsAlert("deamessage"),
                           bsCollapse(id = "collapsedea", open = "DEA plots",
                                      bsCollapsePanel("Genes info", dataTableOutput('deaSE'), style = "default"))),
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
                           box(title = "DEA analysis",width = NULL,
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
                           box(title = "Pathway graphs",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               shinyFilesButton('pathewayexpfile', 'DEA result', 'Please select expression result object',
                                                multiple = FALSE),
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
                           )
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
                           box(title = "Starburst plot",width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = FALSE,
                               box(title = "Data",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                                   bsTooltip("starburstmetfile", "Result from Differential DNA Methylation analysis",
                                             "left"),
                                   shinyFilesButton('starburstmetfile', 'DMR result',
                                                    'Please select DMR_result file',
                                                    multiple = FALSE),
                                   tags$br(),
                                   tags$br(),
                                   bsTooltip("starburstexpfile", "Result from Differential Expression Analysis",
                                             "left"),
                                   shinyFilesButton('starburstexpfile', 'DEA result', 'Please select expression result object',
                                                    multiple = FALSE)),
                               box(title = "Threshold control",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   numericInput("starburstexpFC", "Log FC threshold",
                                                min = 0, max = 10, value = 0, step = 0.05),
                                   numericInput("starburstexFDR", "Expression FDR cut-off",
                                                min = 0, max = 1, value = 0.05, step = 0.001),
                                   numericInput("starburstmetdiff", "Mean DNA methylation difference threshold",
                                                min = 0, max = 1, value = 0, step = 0.05),
                                   numericInput("starburstmetFDR", "Methylation FDR cut-off",
                                                min = 0, max = 1, value = 0.05, step = 0.001)),
                               box(title = "Highlight control",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   checkboxInput("starburstNames", "Show genes names?", value = FALSE, width = NULL),
                                   checkboxInput("starburstNamesFill", "Boxed names?", value = TRUE, width = NULL),
                                   bsTooltip("starburstCircle", "Circle candidate biologically significant genes",
                                             "left"),
                                   checkboxInput("starburstCircle", "Circle genes?", value = TRUE, width = NULL)),

                               box(title = "Colors control",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   colourInput("sbcolInsignigicant", "Insignificant", value = "black"),
                                   colourInput("sbcolUpHypo", "Upregulated & Hypomethylated", value =  "#E69F00"),
                                   colourInput("sbcolDownHypo", "Downregulated & Hypomethylated", value = "#56B4E9"),
                                   colourInput("sbcolHypo", "Hypomethylated", value = "#009E73"),
                                   colourInput("sbcolHyper", "Hypermethylated", value = "red"),
                                   colourInput("sbcolUp", "Upregulated", value =  "#0072B2"),
                                   colourInput("sbcolDown", "Downregulated", value = "#D55E00"),
                                   colourInput("sbcolUpHyper", "Upregulated & Hypermethylated", value = "#CC79A7"),
                                   colourInput("sbcolDownHyper", "Downregulated & Hypermethylated", value = "purple")),
                               box(title = "Size control",width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   sliderInput("starburstwidth", "Plot Width (%)", min = 0, max = 100, value = 100),
                                   sliderInput("starburstheight", "Plot Height (px)", min = 0, max = 1200, value = 800)),
                               bsTooltip("starburstSave", "Save results of significant genes into a csv file",
                                         "left"),
                               checkboxInput("starburstSave", "Save result?", value = FALSE, width = NULL),
                               actionButton("starburstPlot",
                                            "starburst plot",
                                            style = "background-color: #000080;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                            icon = icon("picture-o")))
                    )
                )
        ),
        tabItem(tabName = "elmer",
                fluidRow(
                    column(10,  bsAlert("elmermessage"),
                           bsCollapse(id = "collapelmer", open = "Plots",
                                      bsCollapsePanel("Results table", dataTableOutput('elmerResult'), style = "default"),
                                      bsCollapsePanel("Plots", uiOutput("elmerPlot"), style = "default"))),
                    column(2,
                           box(title = "Create mee", width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE,collapsed = TRUE,
                               shinyFilesButton('elmermetfile', 'Select DNA methylation object', 'Please select DNA methylation object',
                                                multiple = FALSE),
                               tags$br(),
                               tags$br(),
                               shinyFilesButton('elmerexpfile', 'Select expression object', 'Please select gene expression object',
                                                multiple = FALSE),
                               tags$hr(),
                               selectizeInput('elmermeetype',
                                              "Group column",
                                              choices=NULL,
                                              multiple = FALSE),
                               selectizeInput('elmermeesubtype',
                                              "Experiment group",
                                              choices=NULL,
                                              multiple = FALSE),
                               selectizeInput('elmermeesubtype2',
                                              "Control group",
                                              choices=NULL,
                                              multiple = FALSE),
                               tags$hr(),
                               bsTooltip("elmermetnacut", " By default, for the DNA methylation data
                                                            will remove probes with NA values in more than 20% samples and
                                                            remove the anottation data.",
                                         "left"),
                               numericInput("elmermetnacut", "DNA methylation: Cut-off NA samples (%)",
                                            min = 0, max = 1, value = 0.2, step = 0.1),
                               actionButton("elmerpreparemee",
                                            "Create mee object",
                                            style = "background-color: #000080;
                                         color: #FFFFFF;
                                         margin-left: auto;
                                         margin-right: auto;
                                         width: 100%",
                                            icon = icon("floppy-o"))),
                           box(title = "Data", width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
                               shinyFilesButton('elmermeefile', 'Select mee', 'Please select mee object',
                                                multiple = FALSE),
                               tags$br(),
                               tags$br(),
                               shinyFilesButton('elmerresultsfile', 'Select results', 'Please select results object',
                                                multiple = FALSE)),
                           box(title = "Analysis", width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               box(title = "Differently methylated probes", width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   numericInput("elmermetdiff", " DNA methylation difference cutoff",
                                                min = 0, max = 1, value = 0.3, step = 0.05),
                                   numericInput("elmermetpercentage", " percentage",
                                                min = 0, max = 1, value = 0.2, step = 0.01),
                                   numericInput("elmermetpvalue", " pvalue",
                                                min = 0, max = 1, value = 0.01, step = 0.01)
                               ),
                               box(title = "Predict enhancer-gene linkages", width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   numericInput("elmergetpairNumGenes", " Nearby genes",
                                                min = 1, max = 100, value = 20, step = 1),
                                   numericInput("elmergetpairpercentage", " percentage",
                                                min = 0, max = 1, value = 0.2, step = 0.01),
                                   numericInput("elmergetpairpermu", " Number of permuation",
                                                min = 0, max = 10000, value = 10000, step = 1000),
                                   numericInput("elmergetpairpvalue", "Pvalue",
                                                min = 0, max = 1, value = 0.01, step = 0.01),
                                   numericInput("elmergetpairportion", "Portion",
                                                min = 0, max = 1, value = 0.3, step = 0.1),
                                   checkboxInput("elmergetpairdiffExp", "Apply t-test", value = FALSE, width = NULL)
                               ),
                               box(title = "Get enriched motif", width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   numericInput("elmergetenrichedmotifMinIncidence", "Minimum incidence",
                                                min = 0, max = 100, value = 10, step = 1),
                                   numericInput("elmergetenrichedmotifLoweOR", "Lower boundary",
                                                min = 0, max = 5, value = 1.1, step = 0.1)
                               ),
                               box(title = "Identify regulatory TFs", width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   numericInput("elmergetTFpercentage", "Percentage",
                                                min = 0, max = 1, value = 0.2, step = 0.01)
                               ),
                               box(title = "Other options", width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   sliderInput("elmercores", "Cores", step=1,
                                               min = 1, max = parallel::detectCores(), value = 1),
                                   shinyDirButton('elmerFolder', 'Save folder', 'Please select a folder',
                                                  class='shinyDirectories btn-default')
                               ),
                               actionButton("elmerAnalysisBt",
                                            "Run analysis",
                                            style = "background-color: #000080;
                                        color: #FFFFFF;
                                        margin-left: auto;
                                        margin-right: auto;
                                        width: 100%",
                                            icon = icon("flask"))),
                           box(title = "Plots", width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               box(title = "Type selection", width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                                   selectizeInput("elmerPlotType", "Type of plot:",
                                                  c("Scatter plots"="scatter.plot",
                                                    "Schematic Plot"="schematic.plot",
                                                    "Motif enrichment plot"="motif.enrichment.plot",
                                                    "TF ranking plot"="ranking.plot"),
                                                  multiple = FALSE),
                                   #sliderInput("motif.enrichment.plot.or", "OR", min = 1, max = 20, value = 20),
                                   #sliderInput("motif.enrichment.plot.lowerOR", "OR", min = 1, max = 20, value = 1.1),
                                   #sliderInput("motif.enrichment.plot.nbprobes", "NumOfProbes", min = 1, max = 50, value = 10),
                                   selectizeInput("ranking.plot.motif", "Enriched motif:",
                                                  NULL,
                                                  multiple = FALSE),
                                   selectizeInput("ranking.plot.tf", "Highlight TF:",
                                                  NULL,
                                                  multiple = TRUE),
                                   radioButtons("schematic.plot.type", "Schematic plot by:",
                                                c("Genes"="genes",
                                                  "Probes"="probes")),
                                   selectizeInput("schematic.plot.genes", "Genes:",
                                                  NULL,
                                                  multiple = FALSE),
                                   selectizeInput("schematic.plot.probes", "Probes:",
                                                  NULL,
                                                  multiple = FALSE),
                                   radioButtons("scatter.plot.type", "Scatter plot by:",
                                                c("By TF"="tf",
                                                  "By Probes"="probes",
                                                  "By Pair"="pair")),
                                   selectizeInput("scatter.plot.genes", "Genes:",
                                                  NULL,
                                                  multiple = FALSE),
                                   sliderInput("scatter.plot.nb.genes", "Number of genes", min = 1, max = 20, value = 20),
                                   selectizeInput("scatter.plot.probes", "Probes:",
                                                  NULL,
                                                  multiple = FALSE),
                                   selectizeInput("scatter.plot.tf", "TF:",
                                                  NULL,
                                                  multiple = TRUE),

                                   selectizeInput("scatter.plot.motif", "Sites with motif:",
                                                  NULL,
                                                  multiple = FALSE)),
                               box(title = "Size control", width = NULL,
                                   solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                   sliderInput("elmerwidth", "Plot Width (%)", min = 0, max = 100, value = 100),
                                   sliderInput("elmerheight", "Plot Height (px)", min = 0, max = 1200, value = 800)),
                               actionButton("elmerPlotBt",
                                            "Plot",
                                            style = "background-color: #000080;
                                                     color: #FFFFFF;
                                                     margin-left: auto;
                                                     margin-right: auto;
                                                     width: 100%",
                                            icon = icon("picture-o"))),
                           box(title = "Results table", width = NULL,
                               status = "danger",
                               solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,
                               selectizeInput("elmerTableType", "Table to show:",
                                              c("TF"="tf",
                                                "Enriched motifs"="motif",
                                                "Pair probe/gene"="pair",
                                                "Signigicant probes"="sigprobes"
                                              ),
                                              multiple = FALSE))
                    )
                )
        ),
        tabItem(tabName = "config",
                fluidRow(
                    column(1),
                    column(8,
                    box(title = "Configuration options", width = NULL,
                        status = "danger",
                        solidHeader = FALSE, collapsible = FALSE,
                    column(3,
                    bsTooltip("workingDir", "Select where to save outputs","left"),
                    shinyDirButton('workingDir', 'Set working directory', 'Set working directory',
                                   class='shinyDirectories btn-default')),
                    column(9, verbatimTextOutput("wd"))
                    ))))
    )
)

# @title  Client side
# @description Client side - Download data from roadmap project
# @keywords internal
TCGAbiolinksGUIUI <- dashboardPage(
    skin = "blue",
    header,
    sidebar,
    body
)
