library(shiny)
library(shinyFiles)
header <- dashboardHeader(
    title = "biOMICs"
)

sidebar <-  dashboardSidebar(
    sidebarMenu(
        menuItem("Ontology search" , tabName = "ontology", icon = icon("database"))
    )
)

body <-  dashboardBody(
    tagList(
        singleton(tags$head(tags$script(
            paste0('//cdnjs.cloudflare.com/ajax/libs/datatables/',
                   '1.10.5/js/jquery.dataTables.min.js'),
            type = 'text/javascript'))),
        singleton(tags$head(tags$link(href = paste0('//cdn.datatables.net/',
                                                    'tabletools/2.2.3/css/',
                                                    'dataTables.tableTools.css'),
                                      rel = 'stylesheet',type = 'text/css'))),
        singleton(tags$head(tags$link(href = paste0('//cdn.datatables.net',
                                                    '/1.10.5/css/jquery.',
                                                    'dataTables.css'),
                                      rel = 'stylesheet',type = 'text/css'))),
        singleton(tags$head(tags$script(src = paste0('//cdn.datatables.net/',
                                                     'tabletools/2.2.3/js/dat',
                                                     'aTables.tableTools.min.js'),
                                        type = 'text/javascript'))),
        singleton(tags$head(tags$link(rel = "stylesheet", type = "text/css",
                                      href = "biOMICS.css")))
    ),
    tabItems(

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
                               shinyDirButton('directory', 'Folder select', 'Please select a folder',
                                              class='btn action-button', buttonType='warning'),
                               verbatimTextOutput("directorypath"),
                               actionButton("ontSearchBt",
                                            "search",
                                            style = "background-color: #F39C12;
                                            color: #FFFFFF;
                                            margin-left: auto;
                                            margin-right: auto;
                                            width: 100%",
                                            icon = icon("search")),
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
