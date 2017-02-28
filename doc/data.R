## ----message=FALSE, warning=FALSE, include=FALSE-------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE-----------
# Get XML files and parse them
clin.query <- GDCquery(project = "TCGA-READ", data.category = "Clinical", barcode = "TCGA-F5-6702")
GDCdownload(clin.query)
clinical.patient <- GDCprepare_clinic(clin.query, clinical.info = "patient")
clinical.patient.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up")

# Get indexed data
clinical.index <- GDCquery_clinic("TCGA-READ")

## ----echo = TRUE, message = FALSE, warning = FALSE-----------------------
select(clinical.patient,vital_status,days_to_death,days_to_last_followup) %>% datatable
select(clinical.patient.followup, vital_status,days_to_death,days_to_last_followup) %>% datatable
# Vital status should be the same in the follow up table 
filter(clinical.index,submitter_id == "TCGA-F5-6702") %>% select(vital_status,days_to_death,days_to_last_follow_up) %>% datatable

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE-----------
# Get XML files and parse them
recurrent.samples <- GDCquery(project = "TCGA-LIHC",
                             data.category = "Transcriptome Profiling",
                             data.type = "Gene Expression Quantification", 
                             workflow.type = "HTSeq - Counts",
                             sample.type = 	"Recurrent Solid Tumor")$results[[1]] %>% select(cases)
recurrent.patients <- substr(recurrent.samples$cases,1,12)
clin.query <- GDCquery(project = "TCGA-LIHC", data.category = "Clinical", barcode = recurrent.patients)
GDCdownload(clin.query)
clinical.patient <- GDCprepare_clinic(clin.query, clinical.info = "patient") 

## ----echo = TRUE, message = FALSE, warning = FALSE-----------------------
# Get indexed data
GDCquery_clinic("TCGA-LIHC") %>% filter(submitter_id %in% recurrent.patients) %>% 
    select(progression_or_recurrence,days_to_recurrence,tumor_grade) %>% datatable

# XML data
clinical.patient %>% select(bcr_patient_barcode,neoplasm_histologic_grade) %>% datatable


