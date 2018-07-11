# A notification ID
observe({
    shinyDirChoose(input,
                   'IDATfolder',
                   roots = get.volumes(input$workingDir),
                   session = session,
                   restrictions = system.file(package='base'))
})


idat <-  reactive({
    inFile <- input$IDATfolder
    if(class(inFile) != "list") return(NULL)
    baseDir  <- parseDirPath(get.volumes(input$workingDir), input$IDATfolder)
    withProgress(message = 'Reading IDAT files',
                 detail = 'This may take a while...', value = 0, {
                     idat <- tryCatch({
                         read.metharray.exp(base = baseDir,
                                            targets = NULL,
                                            force = TRUE,
                                            recursive = T,
                                            verbose = T )
                     }, warning = function(w){
                         print(w)
                         createAlert(session, "idatAlert", "idatmessage",
                                     title = "Error",
                                     style =  "danger",
                                     content = paste(w), append = FALSE)
                         return(NULL)
                     },error = function(e){
                         print(e)
                         createAlert(session, "idatAlert", "idatmessage",
                                     title = "Error",
                                     style =  "danger",
                                     content = paste(e), append = FALSE)
                         return(NULL)
                     })
                     incProgress(1, detail = "Completed")
                 })
    return(idat)
})
observeEvent(input$IDATfolder, {
    closeAlert(session,"idatmessage")
    output$idattbl <- DT::renderDataTable({
        idat <- idat()
        if (is.null(idat)) return(NULL)
        df <- data.frame(samples = colnames(idat), t(minfi::annotation(idat)))
        createTable(df)
    })
})

observeEvent(input$idatqc, {
    closeAlert(session,"idatmessage")
    idat <- idat()

    withProgress(message = 'Creating ShinyMethyl Summary object',
                 detail = 'This may take a while...', value = 0, {
                     tryCatch({
                         fname <- isolate({input$idatfilename})
                         fname <- gsub(".rda","_ShinyMethylSummary.rda",fname)
                         summary <- shinyMethyl::shinySummarize(idat)
                         incProgress(1, message = "Saving",detail = paste0("As: ", fname))
                         save(summary,file = fname)
                         getPath <- parseDirPath(get.volumes(isolate({input$workingDir})), isolate({input$workingDir}))
                         if (length(getPath) == 0) getPath <- paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI")
                         filename <- file.path(getPath,fname)
                         createAlert(session, "idatAlert", "idatmessage", title = "ShinyMethyl summary data saved", style =  "success",
                                     content =  paste0("Saved in: ", "<br><ul>", paste(filename, collapse = "</ul><ul>"),"</ul>",
                                                       "<br>To visualize QC plot, please use this R command:<ul>load('",filename,"'); shinyMethyl::runShinyMethyl(summary)</ul>"), append = FALSE)
                     }, error = function(e){
                         createAlert(session, "idatAlert", "idatmessage", title = "Error", style =  "error",
                                     content =  paste0("Erro: ", "<br><ul>", paste(e, collapse = "</ul><ul>"),"</ul>"), append = FALSE)
                     })
                 })
})
observeEvent(input$idatnormalize, {
    closeAlert(session,"idatmessage")
    idat <- idat()
    if (is.null(idat)) return(NULL)
    annotation <- as.data.frame(t(minfi::annotation(idat)))

    withProgress(message = 'Normalizing IDAT files',
                 detail = 'This may take a while...', value = 0, {

                     # EPIC has different versions. Later version removed probes we will remove then
                     if(annotation$array == "IlluminaHumanMethylationEPIC") {
                         data("probes2rm")
                         incProgress(0.2, detail = "Filter probes based on B4 file from Illumina")
                         idat <- idat[!rownames(idat) %in% probes2rm,]
                     }

                     incProgress(0.2, detail = "Calculate p-values")
                     ### Step 1.3: Calculate p-values
                     detP <- detectionP(idat, type = "m+u") #failed positions reporting background signal levels

                     ### Step 2: Preprocess:
                     incProgress(0.2, detail = "Noob preprocess")
                     proc <- tryCatch({
                         preprocessNoob(idat, offset = 0,
                                        dyeCorr = TRUE,
                                        verbose = TRUE,
                                        dyeMethod="reference")
                     }, warning = function(w){
                         print(w)
                         createAlert(session, "idatAlert", "idatmessage",
                                     title = "Error",
                                     style =  "danger",
                                     content = paste(w), append = FALSE)
                         return(NULL)
                     },error = function(e){
                         print(e)
                         createAlert(session, "idatAlert", "idatmessage",
                                     title = "Error",
                                     style =  "danger",
                                     content = paste(e), append = FALSE)
                         return(NULL)
                     })
                     if(is.null(proc)) {
                         createAlert(session, "idatAlert", "idatmessage",
                                     title = "Error",
                                     style =  "danger",
                                     content = "Error in preprocessNoob", append = FALSE)
                         return(NULL)
                     }
                     ### Step 3: mask probes that failed p-value detection
                     incProgress(0.2, detail = "Masking probes with detection p-value > 0.05 ")
                     proc.r <- ratioConvert(proc)
                     is.na(assays(proc.r)$Beta) <- (detP[rownames(proc.r), colnames(proc.r)] > 0.05)

                     incProgress(0.2, detail = "Mask probes as recommended by Zhou et al. 2016")
                     # EPIC has different versions. Later version removed probes we will remove then
                     platform <- ifelse(annotation$array == "IlluminaHumanMethylation450k","450K",
                                        ifelse(annotation$array == "IlluminaHumanMethylationEPIC","EPIC","27K"))
                     genome <- ifelse(grepl("hg19",annotation$annotation),"hg19","hg38")

                     met <- ELMER:::getInfiniumAnnotation(genome = genome,plat  = platform)
                     mask <- met[met$MASK_general,]
                     # Remove masked probes, besed on the annotation
                     is.na(assays(proc.r)$Beta) <-   is.na(assays(proc.r)$Beta) | rownames(proc.r) %in% names(mask)

                     ### Step 4: Get beta-values:
                     beta <- as.data.frame(assays(proc.r)$Beta)
                     fname <- isolate({input$idatfilename})
                     getPath <- parseDirPath(get.volumes(isolate({input$workingDir})), isolate({input$workingDir}))
                     if (length(getPath) == 0) getPath <- paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI")
                     fname <- file.path(getPath,fname)

                     incProgress(1, message = "Saving",detail = paste0("As: ", fname))
                     if(tools::file_ext(fname) == "csv"){
                         write.csv(beta,file = fname,row.names = T)
                     } else  if(tolower(tools::file_ext(fname))  %in% c("rda","rdata")){
                         beta <- DataFrame(beta)
                         colnames(beta) <- colnames(idat)
                         beta <- beta[rownames(beta) %in% names(met),,drop=F]
                         rowRanges <- met[rownames(beta)]
                         colData <- DataFrame(Sample=colnames(idat))
                         colData <- cbind(colData,annotation)
                         met <- SummarizedExperiment(assays = SimpleList(beta),rowRanges = rowRanges, colData=colData)
                         save(met,file = fname)
                     } else {
                         createAlert(session,
                                     "idatAlert",
                                     "idatmessage",
                                     title = "Error in saving",
                                     style =  "danger",
                                     content = "File name has to be .csv or .rda. Saved as Idat.rda", append = FALSE)
                         save(beta,file = file.path(getPath,"Idat.rda"))
                     }

                     createAlert(session, "idatAlert", "idatmessage", title = "Processed data saved", style =  "success",
                                 content =  paste0("Saved in: ", "<br><ul>", paste(fname, collapse = "</ul><ul>"),"</ul>"), append = FALSE)
                 })

})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Download button
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observe({
    output$idatdownloadDataBut <- downloadHandler(
        filename <-   function() {
            as.character( isolate({input$idatfilename}))
        },
        content <- function(file){
            file.copy( isolate({input$idatfilename}),file)
        }
    )
})
