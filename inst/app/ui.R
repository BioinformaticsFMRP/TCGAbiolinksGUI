.header <- dashboardHeader(
  title = "biOMICs"
)

.sidebar <-  dashboardSidebar(
  sidebarMenu(
    menuItem ("Encode" , tabName = "encode" , icon = icon("download")),
    menuItem ("Roadmap", tabName = "roadmap", icon = icon("download")),
    menuItem ("Roadmap table", tabName = "roadmap_table", icon = icon("download")),
    menuItem ("TCGA"   , tabName = "tcga"   , icon = icon("database"), badgeLabel = "new", badgeColor = "green")
  )
)

.body <-  dashboardBody(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "biOMICS.css")
  ),
  tabItems(
    # First tab content
    tabItem(tabName = "roadmap",

            fluidRow(
              column(8,
                     box(title = "Advanced search",width = NULL, status = "warning",
                         solidHeader = TRUE, collapsible = FALSE,
                         textInput("rmapSearch", label = "Search Term", value = ""),
                         #textInput("rmapProject", label = "Project", value = "roadmap epigenomics"),
                         #textInput("rmapType", label = "Type", value = "gsm"),
                         selectInput('rmapExp', 'Experiments filter',
                                     roadmap$exp,
                                     multiple = TRUE, selectize = TRUE),
                         selectInput('rmapSamples', 'Samples filter',
                                     roadmap$samples,
                                     multiple = TRUE, selectize = TRUE),
                         actionButton("rmapSearchBt",
                                      "search",
                                      style = "background-color: #F39C12;color: #FFFFFF;
                                      margin-left: auto;margin-right: auto;width: 49%",
                                      icon = icon("search")),
                         actionButton("rmapDownloadBt",
                                      "Download",
                                      style = "background-color: #F39C12;color: #FFFFFF;
                                      margin-left: auto;margin-right: auto;width: 49%",
                                      icon = icon("download"))
                     )

              ),
              column(4,
                     valueBoxOutput("savedPath", width = NULL),
                     uiOutput("statusBox"),
                     uiOutput("savedFiles"))
            )
    ),
    tabItem(tabName = "roadmap_table",

            fluidRow(
              column(1),
              column(10, dataTableOutput('rmaptbl')),
              column(1)
            ),
            fluidRow(
              column(4),
              column(4,
              actionButton("rmapTableDownloadBt",
                           "Download selected rows",
                           style = "background-color: #F39C12;color: #FFFFFF;
                                      margin-left: auto;margin-right: auto;width: 100%",
                           icon = icon("download")),
              textOutput('rmapTableLink')
              ),
              column(4)

            )
    ),
    # Second tab content
    tabItem(tabName = "encode",
            fluidRow(
              column(8,
                     box(title = "Advanced search",width = NULL, status = "warning",
                         solidHeader = TRUE, collapsible = TRUE,
                         textInput("encodeSearch", label = "Search Term", value = ""),
                         selectInput('encodeAssay', 'Assay filter',
                                     list("ChIP-seq",
                                          "RNA-seq",
                                          "RRBS",
                                          "RIP-chip",
                                          "DNase-seq",
                                          "Repli-chip"),
                                     multiple = TRUE, selectize = TRUE),
                         # @param iTarget - "transcription factor" "tag" ...
                         selectInput("encodeTarget",
                                     label =  "Target filter",
                                     choices = list("transcription factor",
                                                    "RNA binding protein",
                                                    "control",
                                                    "histone modification",
                                                    "tag"),
                                     multiple = TRUE, selectize = TRUE),
                         # @param iFType  - "bam" "bigWig" "bed_broadPeak" ...
                         selectInput("encodeFtype",
                                     label = "File Type filter",
                                     choices = list("bam",
                                                    "bigWig",
                                                    "bed_broadPeak",
                                                    "broadPeak",
                                                    "narrowPeak",
                                                    "bed_narrowPeak",
                                                    "bed",
                                                    "bigBed",
                                                    "fastq"),
                                     multiple = TRUE, selectize = TRUE),
                         # @param iSample - "tissue" "primary cell"
                         selectInput("encodeSample",
                                     label = "Sample filter",
                                     choices = list("tissue",
                                                    "stem cell",
                                                    "immortalized cell line",
                                                    "in vitro differentiated cells",
                                                    "induced pluripotent stem cell line",
                                                    "primary cell"),
                                     multiple = TRUE, selectize = TRUE),
                         # @param assembly - "hg19" "mm9"
                         selectInput("encodeAssembly",
                                     label ="Assembly Filter",
                                     choices = list("hg19",
                                                    "dm3",
                                                    "mm9"),
                                     multiple = TRUE, selectize = TRUE),
                         actionButton("encodeDownloadBt",
                                      "Download",
                                      style = "background-color: #F39C12;color: #FFFFFF;
                                      margin-left: auto;margin-right: auto;width: 100%",
                                      icon = icon("download"))
                     )
              ),
              column(4,
                     valueBoxOutput("savedPath2", width = NULL)
              )
            ),
            fluidRow(
              column(1),
              column(10, dataTableOutput('tblEncode')),
              column(1)
            )
    ),
    tabItem(tabName = "tcga",
            fluidRow(
              column(1),
              column(10,
                     box(title = "Get bar code", width = NULL, status = "warning",
                         solidHeader = TRUE, collapsible = FALSE,
                         fileInput('file1', 'Upload TCGA file',
                                   accept = c(
                                     'text/csv',
                                     'text/comma-separated-values',
                                     '.csv'
                                   )
                         ),
                         fileInput('file2', 'Upload barcode filter',
                                   accept = c(
                                     'text/plain'
                                   )),
                         downloadButton("getTcgaBarCode",
                                      "Download",
                                      class = "btn-block btn-warning"
                                      )

                     )),
              column(1)
            )
    )
  )
)

# @title  Client side
# @description Client side - Download data from roadmap project
# @keywords internal
.biOMICsUI <- dashboardPage(
  skin = "blue",
  .header,
  .sidebar,
  .body
)
