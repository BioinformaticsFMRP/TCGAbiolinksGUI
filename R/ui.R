.header <- dashboardHeader(
  title = "biOMICs",
  dropdownMenu(type = "messages",
               messageItem(
                 from = "Sales Dept",
                 message = "Sales are steady this month."
               ),
               messageItem(
                 from = "New User",
                 message = "How do I register?",
                 icon = icon("question"),
                 time = "13:45"
               ),
               messageItem(
                 from = "Support",
                 message = "The new server is ready.",
                 icon = icon("life-ring"),
                 time = "2014-12-01"
               )
  ),
  dropdownMenu(type = "notifications",
               notificationItem(
                 text = "5 new users today",
                 icon("users")
               ),
               notificationItem(
                 text = "12 items delivered",
                 icon("truck"),
                 status = "success"
               ),
               notificationItem(
                 text = "Server load at 86%",
                 icon = icon("exclamation-triangle"),
                 status = "warning"
               )
  ),

  dropdownMenu(type = "tasks", badgeStatus = "success",
               taskItem(value = 40, color = "green",
                        "Documentation"
               ),
               taskItem(value = 70, color = "aqua",
                        "Project roadmap"
               ),
               taskItem(value = 30, color = "yellow",
                        "UI roadmap"
               ),
               taskItem(value = 5, color = "red",
                        "Overall project"
               )
  )

)

.sidebar <-  dashboardSidebar(
  sidebarMenu(
    menuItem("Encode", tabName = "encode", icon = icon("download")),
    menuItem("Roadmap", tabName = "roadmap", icon = icon("download"),badgeLabel = "new", badgeColor = "green")

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
                         textInput("rmapProject", label = "Project", value = "roadmap epigenomics"),
                         textInput("rmapType", label = "Type", value = "gsm"),
                         actionButton("rmapDownloadBt",
                                      "Download",
                                      style = "background-color: #F39C12;color: #FFFFFF;
                                      margin-left: auto;margin-right: auto;width: 100%",
                                      icon = icon("download"))
                     )
              ),
              column(4,
                     valueBoxOutput("savedPath", width = NULL),
                     uiOutput("statusBox"),
                     uiOutput("savedFiles"))
            ),
            fluidRow(
              column(1),
              column(10, DT::dataTableOutput('tbl')),
              column(1)
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
              column(10, DT::dataTableOutput('tblEncode')),
              column(1)
            )
    )
  )
)

#' @title  Client side
#' @description Client side - Download data from roadmap project
#' @name biOMICsUI
#' @keywords internal
#'
biOMICsUI <- dashboardPage(
  skin = "blue",
  .header,
  .sidebar,
  .body
)
