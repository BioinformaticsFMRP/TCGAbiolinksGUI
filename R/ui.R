library("shiny")

shinyUI(fluidPage(
  titlePanel('Downloading ENCODE Data'),
  sidebarLayout(
    sidebarPanel(
      # Copy the line below to make a text input box
      textInput("text", label = h3("Search Term"), value = "Enter text..."),
  

      hr(),
     
      #' @param assay - "ChIP-seq" "RIP-chip" "Repli-chip" 
    checkboxGroupInput("assay", label = h3("Assay"), 
                       choices = list("ChIP-seq", "RIP-chip", "Repli-chip"),
                       selected = NULL),
    
    hr(),
    fluidRow(column(3, verbatimTextOutput("assay"))),
    
    hr(),
 
    #' @param iTarget - "transcription factor" "RNA binding protein" "tag" "histone modification"
    checkboxGroupInput("target", label = h3("Target"), 
                       choices = list("transcription factor",  "RNA binding protein",  "histone modification","tag"),
                       selected = NULL),
    
    hr(),
    fluidRow(column(3, verbatimTextOutput("target"))),
    
    
    #' @param iFType  - "bam" "bigWig" "bed_broadPeak" "broadPeak" "fastq"
    checkboxGroupInput("ftype", label = h3("File Type"), 
                       choices = list("bam","bigWig","bed_broadPeak","broadPeak","fastq"),
                       selected = NULL),
    
    hr(),
    fluidRow(column(3, verbatimTextOutput("ftype"))),
    
    #' @param iSample - "tissue" "primary cell"
    checkboxGroupInput("sample", label = h3("Sample"), 
                       choices = list("tissue","primary cell"),
                       selected = NULL),
    
    hr(),
    fluidRow(column(3, verbatimTextOutput("ftype"))),
    #' @param assembly - "hg19" "mm9"
    
    
    
    hr(),
    downloadButton('downloadData', 'Download')
     
    ),
    mainPanel(
      tableOutput('table')
    )
  )
))