#' @title  Client side 
#' @description Client side - Download data from ENCODDE project
#' @name encodeDownloaderUI
#' @keywords internal

shinyUI <- fluidPage(
  
  tags$head(tags$style(HTML("
    .shiny-text-output {
      background-color:#fff;
    }
  "))),
  
  h1("Downloading", span("ENCODE Data", style = "font-weight: 300"), 
     style = "font-family: 'Source Sans Pro';
        color: #fff; text-align: center;
        background-image: url('texturebg.png');
        background-color: #4187C5; 
        padding: 20px"),
  br(),
  
  
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
