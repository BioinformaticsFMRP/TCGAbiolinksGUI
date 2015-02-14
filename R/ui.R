header <- dashboardHeader(
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
               taskItem(value = 90, color = "green",
                        "Documentation"
               ),
               taskItem(value = 17, color = "aqua",
                        "Project X"
               ),
               taskItem(value = 75, color = "yellow",
                        "Server deployment"
               ),
               taskItem(value = 80, color = "red",
                        "Overall project"
               )
  )

  )

sidebar <-  dashboardSidebar(
  sidebarMenu(
    menuItem("Encode", tabName = "encode", icon = icon("download")),
    menuItem("Roadmap", tabName = "roadmap", icon = icon("download"),badgeLabel = "new", badgeColor = "green")

    )
)

body <-  dashboardBody(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "biOMICS.css")
  ),
  tabItems(
    # First tab content
    tabItem(tabName = "roadmap",

            fluidRow(
              column(4,
                     box(title = "Advanced search",width = NULL, status = "warning",
                         solidHeader = TRUE, collapsible = TRUE,
                         textInput("search", label = "Search Term", value = ""),
                         textInput("project", label ="Project", value = "roadmap epigenomics"),
                         textInput("type", label = "Type", value = "gsm")
                     )
              ),
              column(6,

                     box(
                       title = "Results", width = NULL, background = "light-blue",
                       verbatimTextOutput("nbresult")
                     )
              )
            ),

            fluidRow(
              column(12,
                     wellPanel(align = "center",
                               style = "background-color: #4187C5",
                               actionButton("downloadBt",
                                            "Download",
                                            icon = icon("download")
                               )
                     )
              )
            )
    ),

    # Second tab content
    tabItem(tabName = "encode",
            fluidRow(

              column(4,
                     wellPanel(
                       textInput("search", label = h3("Search Term"), value = "")
                     )
              ),

              column(4,
                     wellPanel(

                       # @param assay - "ChIP-seq" "RIP-chip" "Repli-chip"
                       checkboxGroupInput("assay",
                                          label = h3("Assay"),
                                          choices = list("ChIP-seq",
                                                         "RNA-seq",
                                                         "RRBS",
                                                         "RIP-chip",
                                                         "DNase-seq",
                                                         "Repli-chip"),
                                          selected = NULL)
                     )
              ),
              column(4,
                     wellPanel(

                       # @param iTarget - "transcription factor" "tag" ...
                       checkboxGroupInput("target",
                                          label = h3("Target"),
                                          choices = list("transcription factor",
                                                         "RNA binding protein",
                                                         "control",
                                                         "histone modification",
                                                         "tag"),
                                          selected = NULL)
                     )
              )
            ),


            fluidRow(
              column(4,
                     wellPanel(
                       # @param iFType  - "bam" "bigWig" "bed_broadPeak" ...
                       checkboxGroupInput("ftype",
                                          label = h3("File Type"),
                                          choices = list("bam",
                                                         "bigWig",
                                                         "bed_broadPeak",
                                                         "broadPeak",
                                                         "narrowPeak",
                                                         "bed_narrowPeak",
                                                         "bed",
                                                         "bigBed",
                                                         "fastq"),
                                          selected = NULL)

                     )
              ),

              column(4,
                     wellPanel(
                       # @param iSample - "tissue" "primary cell"
                       checkboxGroupInput("sample",
                                          label = h3("Sample"),
                                          choices = list("tissue",
                                                         "stem cell",
                                                         "immortalized cell line",
                                                         "in vitro differentiated cells",
                                                         "induced pluripotent stem cell line",
                                                         "primary cell"),
                                          selected = NULL)
                     )
              ),

              column(4,
                     wellPanel(
                       # @param assembly - "hg19" "mm9"
                       checkboxGroupInput("assembly",
                                          label = h3("Assembly"),
                                          choices = list("hg19",
                                                         "dm3",
                                                         "mm9"),
                                          selected = NULL)
                     )
              )
            ),

            fluidRow(
              column(12,
                     wellPanel(align = "center",
                               style = "background-color: #4187C5",
                               actionButton("downloadBt",
                                            "Download",
                                            icon = icon("download")
                               )
                     )
              )
            ),
            fluidRow(column(5, offset = 4, verbatimTextOutput("value")))

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
  header,
  sidebar,
  body
)
