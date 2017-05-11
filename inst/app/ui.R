suppressPackageStartupMessages({
    library(shiny)
    library(shinyFiles)
    library(TCGAbiolinks)
    library(shinyBS)
    library(shinyjs)
    library(SummarizedExperiment)
    library(pathview)
    library(reshape2)
    library(shinydashboard)
    library(clusterProfiler)
})

data(paths.hsa)
pathways.id <- names(paths.hsa)
names(pathways.id) <- unname(paths.hsa)
pathways.id <- pathways.id[sort(unname(paths.hsa))]
menu.icon <- "arrow-circle-right"
ui.path <- ifelse(system.file("app", package = "TCGAbiolinksGUI") == "",
                  "ui",
                  file.path(system.file("app", package = "TCGAbiolinksGUI"),"ui"))

#' busyIndicator
#'
# This is a function to indicate the work is in progress, it was created for the plots
# that rendering were taking long and withprogress was not working.
# @param text The text to show
# @param wait The amount of time to wait before showing the busy indicator. The
#   default is 1000 which is 1 second.
#
# @export
busyIndicator <- function(text = "Working in progress...") {
    div(
        id = 'busyModal', class = 'modal', role = 'dialog', 'data-backdrop' = 'static',
        div(
            class = 'modal-dialog modal-sm',
            div(id = 'modal-content-busy',
                class = 'modal-content',
                div(class = 'modal-header', h4(class = 'modal-title', text)),
                div(class = 'modal-body', p(h2(HTML('<i class="fa fa-cog fa-spin"></i>'))))
            )
        )
    )
}

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
        tags$textarea(id = inputId,
                      class = "inputtextarea",
                      rows = nrows,
                      cols = ncols,
                      as.character(value))
    )
}


header <- dashboardHeader(
    title = "TCGAbiolinksGUI",
    titleWidth = 250
)
header$children[[2]]$children <-  tags$a(href='http://bioconductor.org/packages/TCGAbiolinksGUI/',
                                         tags$img(src='logo_no_text.png',height='50',width='220'))

sidebar <-  dashboardSidebar(
    width = 250,
    sidebarMenu(
        tags$hr(class="lineData"),
        menuItem("Get GDC data", icon = icon("database"),
                 menuSubItem("Molecular data" , tabName = "tcgaSearch", icon = icon("database")),
                 menuSubItem("Mutation data" , tabName = "tcgaMutation", icon = icon("database")),
                 menuSubItem("Clinical data" , tabName = "tcgaClinical", icon = icon("database")),
                 menuSubItem("Subtype data" , tabName = "tcgaSubtype", icon = icon("database"))
        ),
        menuItem("Manage SummarizedExperiment" , tabName = "seedit", icon = icon("pencil")),
        tags$hr(class="lineAnalysis"),
        menuItem("Clinical analysis", icon = icon("flask"),
                 menuSubItem("Survival plot" , tabName = "tcgasurvival", icon = icon("picture-o"))
        ),

        menuItem("Epigenetic analysis", icon = icon("flask"),
                 menuSubItem("Differential methylation analysis" , tabName = "dmr", icon = icon("flask")),
                 menuSubItem("Volcano plot" , tabName = "volcano", icon = icon("picture-o")),
                 menuSubItem("Heatmap plot" , tabName = "heatmap", icon = icon("picture-o")),
                 menuSubItem("Mean DNA methylation plot" , tabName = "meanmet", icon = icon("picture-o"))
        ),
        menuItem("Transcriptomic analysis", icon =  icon("flask"),
                 menuSubItem("Differential expression analysis" , tabName = "dea", icon = icon("flask")),
                 menuSubItem("Volcano plot" , tabName = "volcano", icon = icon("picture-o")),
                 menuSubItem("Heatmap plot" , tabName = "heatmap", icon = icon("picture-o")),
                 menuSubItem("Pathway visualization" , tabName = "pathview", icon = icon("picture-o")),
                 menuSubItem("Enrichment analysis" , tabName = "ea", icon = icon("flask")),
                 menuSubItem("Network inference" , tabName = "netinf", icon = icon("flask"))
        ),
        menuItem("Genomic analysis", icon = icon("flask"),
                 menuSubItem("OncoPrint plot" , tabName = "tcgaOncoPrint", icon = icon("picture-o"))
        ),
        tags$hr(class="lineIntegrative"),
        menuItem("Starburst plot" , tabName = "starburst", icon = icon("picture-o")),
        menuItem("ELMER" , icon = icon("flask"),
                 menuSubItem("Analysis" , tabName = "elmeranalysis", icon = icon("flask")),
                 menuSubItem("Visualize results" , tabName = "elmerresults", icon = icon("picture-o"))
        ),
        tags$hr(class="lineConfig"),
        menuItem("Configuration", tabName = "config", icon = icon("cogs")),
        tags$hr(class="lineDoc"),
        menuItem("Tutorial/Vignettes", icon = icon("book"),
                 menuSubItem("TCGAbiolinksGUI Manual" , href = "http://bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinksGUI/inst/doc/index.html", icon = icon("external-link")),
                 menuSubItem("TCGAbiolinks Manual" , href = "https://www.bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/index.html", icon = icon("external-link")),
                 menuSubItem("ELMER Manual" , href = "https://www.bioconductor.org/packages/3.3/bioc/vignettes/ELMER/inst/doc/vignettes.pdf", icon = icon("external-link"))
        ),
        menuItem("References", icon = icon("file-text-o"), tabName = "references"
                 #menuSubItem("TCGAbiolinks" , href = "https://doi.org/10.1093/nar/gkv1507", icon = icon("external-link")),
                 #menuSubItem("ELMER" , href = "https://doi.org/10.1186/s13059-015-0668-3", icon = icon("external-link"))
        ), div(id = "greetbox-outer",
               menuItem("Loading...", tabName = "welcome", icon = icon("spinner"), selected = T)))
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
                                      href = "TCGAbiolinksGUI.css"))),
        singleton(tags$head(tags$script(src = 'events.js')))
    ),
    tabItems(
        source(file.path(ui.path, "index.R"),  local = TRUE)$value,
        source(file.path(ui.path, "getmolecular.R"),  local = TRUE)$value,
        source(file.path(ui.path, "manageSE.R"),  local = TRUE)$value,
        source(file.path(ui.path, "getclinical.R"),  local = TRUE)$value,
        source(file.path(ui.path, "getmutation.R"),  local = TRUE)$value,
        source(file.path(ui.path, "getsubtype.R"),  local = TRUE)$value,
        source(file.path(ui.path, "oncoprint.R"),  local = TRUE)$value,
        source(file.path(ui.path, "volcano.R"),  local = TRUE)$value,
        source(file.path(ui.path, "heatmap.R"),  local = TRUE)$value,
        source(file.path(ui.path, "dmr.R"),  local = TRUE)$value,
        source(file.path(ui.path, "meanMet.R"),  local = TRUE)$value,
        source(file.path(ui.path, "ea.R"),  local = TRUE)$value,
        source(file.path(ui.path, "survival.R"),  local = TRUE)$value,
        source(file.path(ui.path, "dea.R"),  local = TRUE)$value,
        source(file.path(ui.path, "pathview.R"),  local = TRUE)$value,
        source(file.path(ui.path, "starburst.R"),  local = TRUE)$value,
        source(file.path(ui.path, "elmer_analysis.R"),  local = TRUE)$value,
        source(file.path(ui.path, "elmer_results.R"),  local = TRUE)$value,
        source(file.path(ui.path, "config.R"),  local = TRUE)$value,
        source(file.path(ui.path, "references.R"),  local = TRUE)$value,
        source(file.path(ui.path, "getinference.R"),  local = TRUE)$value

    )
)

# @title  Client side
# @description Client side - Download data from roadmap project
# @keywords internal
shinyUI(
    bootstrapPage(
        useShinyjs(),
        div(id = "loading-content",
            img(src = "progress.gif"),
            tags$br(),
            "Loading TCGAbiolinksGUI..."
        ),
        dashboardPage(
            skin = "blue",
            header,
            sidebar,
            body),
        busyIndicator() # Add rendering in progress...
    )
)
