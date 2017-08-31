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
                                             class='shinyDirectories btn-default'),
                              downloadButton('downloadDataBut', 'Download', class = "butt2"),
                              # making the font italics this time
                              tags$head(tags$style(".butt2{background-color:navy;} .butt2{color: white;} .butt2{width: 100%;} .butt2{margin-top: 11pt;} .butt2{border-radius: 5px;}"))
                       ),
                       column(9, verbatimTextOutput("wd"),
                              shinyFilesButton('downloadfile', 'Select file to download', 'Select file to download',
                                               multiple = FALSE))
                   )
            )
        )
)
