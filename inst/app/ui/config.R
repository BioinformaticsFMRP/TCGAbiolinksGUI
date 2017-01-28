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
                   )))
)
