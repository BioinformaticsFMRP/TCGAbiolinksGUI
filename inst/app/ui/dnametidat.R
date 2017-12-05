tabItem(tabName = "dnametidat",
        fluidRow(
            column(8, bsAlert("idatAlert"),
                   bsCollapse(id = "collapseIdat", open = "Samples identified",
                              bsCollapsePanel("Samples identified", DT::dataTableOutput('idattbl'), style = "default")
                   )),
            column(4,
                   box(title = "IDAT normalization",width = NULL,
                       status = "danger",
                       solidHeader = FALSE, collapsible = FALSE,
                       shinyDirButton('IDATfolder', 'Select idat folder', 'Please select a folder', TRUE),
                       tags$hr(),
                       selectizeInput('idatmetPlatform',
                                      label = 'DNA methylation platform',
                                      choices = c("Illumina Human Methylation 450" = "450K",
                                                  "Illumina Human Methylation 27" = "27K",
                                                  "Infinium MethylationEPIC"= "EPIC"),
                                      multiple = FALSE),
                       bsTooltip("idatmetPlatform", "Used to mask probes as recommended by Zhou et al. 2016 (http://zwdzwd.io/InfiniumAnnotation/current/)",
                                 "right"),
                       selectizeInput('idatgenome',
                                      label = 'Genome of reference',
                                      choices = c("hg38"="hg38","hg19"="hg19"),
                                      multiple = FALSE),
                       bsTooltip("idatgenome", "Used to mask probes as recommended by Zhou et al. 2016 (http://zwdzwd.io/InfiniumAnnotation/current/)",
                                 "right"),
                       textInput("idatfilename", label = "Filename", value = "normalizedIdat.csv"),
                       bsTooltip("idatfilename", "Extensions: csv (flat table) or rda (summarized Experiment)",
                                 "left"),
                       actionButton("idatnormalize",
                                    "Normalize and save",
                                    style = "background-color: #000080;
                     color: #FFFFFF;
                     margin-left: auto;
                     margin-right: auto;
                     width: 100%",
                                    icon = icon("flask")),
                       tags$hr(),
                       downloadButton('idatdownloadDataBut', 'Download normalized results',
                                      style = "background-color: #000080;
                       color: #FFFFFF;
                       margin-left: auto;
                       margin-right: auto;
                       width: 100%")),
                   box(title = "Glioma classification",width = NULL,
                       status = "danger",
                       solidHeader = FALSE, collapsible = FALSE,
                       shinyFilesButton(id = 'classifyObj',
                                        label = 'Select DNA methylation (.rda) file',
                                        title = 'Please select a DNA methylation (.rda) file',
                                        multiple = FALSE),
                       tags$br(),
                       tags$br(),
                       actionButton("idatClassify",
                                    "Classify into gliomas molecular subtypes",
                                    style = "background-color: #000080;
                                    color: #FFFFFF;
                                    margin-left: auto;
                                    margin-right: auto;
                                    width: 100%",
                                    icon = icon("flask"))
                   )))
)
