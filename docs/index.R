## ---- eval = FALSE-------------------------------------------------------
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("TCGAbiolinksGUI", dependencies = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  source("https://bioconductor.org/biocLite.R")
#  deps <- c("pathview","clusterProfiler","ELMER", "DO.db","GO.db", "ComplexHeatmap","EDASeq", "TCGAbiolinks")
#  for(pkg in deps)  if (!pkg %in% installed.packages()) biocLite(pkg, dependencies = TRUE)
#  deps <- c("devtools","shape","shiny","readr","googleVis","shinydashboard","shinyFiles","shinyjs","shinyBS")
#  for(pkg in deps)  if (!pkg %in% installed.packages())  install.packages(pkg,dependencies = TRUE)
#  devtools::install_github("BioinformaticsFMRP/TCGAbiolinksGUI")

## ---- eval = FALSE-------------------------------------------------------
#  library(TCGAbiolinksGUI)
#  TCGAbiolinksGUI()

## ----table2, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "  
| Menu                    | Sub-menu                          | Button             | Data input  |
|---------------------------------|-----------------------------------|------------------------------------|-----------------------------------------------------------------------------------------------------------------------|
| Clinical analysis       | Profile Plot                      | Select file        | A table with at least two categorical columns  |
| Clinical analysis       | Survival Plot                     | Select file        | A table with at least the following columns: days_to_death, days_to_last_followup and one column with a group |
| Epigenetic analysis     | Differential methylation analysis | Select data (.rda) | A summarizedExperiment object |
| Epigenetic analysis     | Volcano Plot                      | Select results     | A csv file with the following pattern: DMR_results_GroupCol_group1_group2_pcut_1e-30_meancut_0.55.csv  (Where GroupCol, group1, group2 are the names of the columns selected in the  DMR steps. |
| Epigenetic analysis     | Heatmap plot                      | Select file        | A summarizedExperiment object  |
| Epigenetic analysis     | Heatmap plot                      | Select results     | Same as Epigenetic analysis >Volcano Plot > Select results |
| Epigenetic analysis     | Mean DNA methylation              | Select file        | A summarizedExperiment object |
| Transcriptomic Analysis | Volcano Plot                      | Select results     | A csv file with the following pattern: DEA_results_GroupCol_group1_group2_pcut_1e-30_meancut_0.55.csv (Where GroupCol, group1, group2 are the names of the columns selected in the DEA steps. |
| Transcriptomic Analysis | Heatmap plot                      | Select file        | A summarizedExperiment object  |
| Transcriptomic Analysis | OncoPrint plot                      | Select MAF file        | A MAF file (columns needed: Hugo_Symbol,Tumor_Sample_Barcode,Variant_Type) |   |
| Transcriptomic Analysis | OncoPrint plot                      | Select Annotation file        | A file with at least the following columns: bcr_patient_barcode  |    |
|  Integrative analysis   | Starburst plot                      | DMR result        | A csv file with the following pattern: DMR_results_GroupCol_group1_group2_pcut_1e-30_meancut_0.55.csv (Where GroupCol, group1, group2 are the names of the columns selected in the DMR steps.  |
|  Integrative analysis   | Starburst plot                      | DEA result        | A csv file with the following pattern: DEA_results_GroupCol_group1_group2_pcut_1e-30_meancut_0.55.csv (Where GroupCol, group1, group2 are the names of the columns selected in the DEA steps.  |
|  Integrative analysis   | ELMER                      | Create mee > Select DNA methylation object         | An rda file with a summarized Experiment object |   
|  Integrative analysis   | ELMER                      | Select results > Select expression object         |  An rda file with the RNAseq data frame |   
|  Integrative analysis   | ELMER                      | Select mee         | An rda file with a mee object |   
|  Integrative analysis   | ELMER                      | Select results         | An rda file with the results of the ELMER analysis |   
"
cat(tabl) 

## ---- include=FALSE------------------------------------------------------
library(TCGAbiolinksGUI)

## ------------------------------------------------------------------------
sessionInfo()

