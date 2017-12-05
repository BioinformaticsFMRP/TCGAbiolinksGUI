# A notification ID
observe({
    shinyDirChoose(input,
                   'IDATfolder',
                   roots = get.volumes(input$workingDir),
                   session = session,
                   restrictions = system.file(package='base'))
})


idat <-  reactive({
    closeAlert(session,"idatAlert")
    inFile <- input$IDATfolder
    if (is.null(inFile)) return(NULL)
    baseDir  <- parseDirPath(get.volumes(input$workingDir), input$IDATfolder)
    withProgress(message = 'Reading IDAT files',
                 detail = 'This may take a while...', value = 0, {
                     idat <- tryCatch({
                         read.metharray.exp(base = baseDir,
                                            targets = NULL,
                                            force = TRUE,
                                            recursive = T,
                                            verbose=T)
                     }, warning = function(w){
                         print(w)
                         createAlert(session, "idatmessage", "idatAlert",
                                     title = "Error",
                                     style =  "danger",
                                     content = paste(w), append = FALSE)
                         return(NULL)
                     },error = function(e){
                         print(e)
                         createAlert(session, "idatmessage", "idatAlert",
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
    output$idattbl <- DT::renderDataTable({
        idat <- idat()
        if (is.null(idat)) return(NULL)
        df <- as.data.frame(colnames(idat))
        colnames(df) <- "IDAT"
        df
    })
})

observeEvent(input$idatnormalize, {
    closeAlert(session,"idatAlert")
    idat <- idat()
    if (is.null(idat)) return(NULL)

    # Quality control report (Visualization):
    # summary <- shinySummarize(idat)
    # runShinyMethyl(summary)

    withProgress(message = 'Normalizing IDAT files',
                 detail = 'This may take a while...', value = 0, {

                     # EPIC has different versions. Later version removed probes we will remove then
                     if(isolate({input$idatmetPlatform}) == "EPIC") {
                         incProgress(0.2, detail = "Filter probes based on B4 file from Illumina")
                         idat <- idat[!rownames(idat) %in% prob_to_remove,]
                     }

                     incProgress(0.2, detail = "Calculate p-values")
                     ### Step 1.3: Calculate p-values
                     detP <- detectionP(idat, type = "m+u") #failed positions reporting background signal levels
                     # table(detP > 0.05)


                     ### Step 2: Preprocess:
                     incProgress(0.2, detail = "Noob preprocess")
                     proc <- tryCatch({
                         preprocessNoob(idat, offset = 0,
                                        dyeCorr = TRUE,
                                        verbose = TRUE,
                                        dyeMethod="reference")
                     }, warning = function(w){
                         print(w)
                         createAlert(session, "idatmessage", "idatAlert",
                                     title = "Error",
                                     style =  "danger",
                                     content = paste(w), append = FALSE)
                         return(NULL)
                     },error = function(e){
                         print(e)
                         createAlert(session, "idatmessage", "idatAlert",
                                     title = "Error",
                                     style =  "danger",
                                     content = paste(e), append = FALSE)
                         return(NULL)
                     })
                     if(is.null(proc)) {
                         createAlert(session, "idatmessage", "idatAlert",
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
                     met <- ELMER:::getInfiniumAnnotation(genome = isolate({input$idatgenome}),plat  = isolate({input$idatmetPlatform}))
                     mask <- met[met$MASK.general,]
                     # Remove masked probes, besed on the annotation
                     is.na(assays(proc.r)$Beta) <-   is.na(assays(proc.r)$Beta) | rownames(proc.r) %in% names(mask)

                     ### Step 4: Get beta-values:
                     beta <- as.data.frame(assays(proc.r)$Beta)
                     fname <- isolate({input$idatfilename})
                     incProgress(1, message = "Saving",detail = paste0("As: ", fname))
                     if(tools::file_ext(fname) == "csv"){
                         write.csv(beta,file = fname,row.names = T)
                     } else  if(tolower(tools::file_ext(fname))  %in% c("rda","rdata")){
                         beta <- DataFrame(beta)
                         beta <- beta[rownames(beta) %in% names(met),]
                         rowRanges <- met[rownames(beta)]
                         colData <- DataFrame(Sample=colnames(beta))
                         met <- SummarizedExperiment(assays = SimpleList(beta),rowRanges = rowRanges, colData=colData)
                         save(met,file = fname)
                     } else {
                         createAlert(session, "idatmessage",
                                     "idatAlert",
                                     title = "Error in saving",
                                     style =  "danger",
                                     content = "File name has to be .csv or .rda. Saved as Idat.rda", append = FALSE)
                         save(beta,file = "Idat.rda")
                     }
                     createAlert(session, "idatmessage", "idatAlert", title = "idatnormalized data saved", style =  "success",
                                 content = "Click in the download button", append = FALSE)
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

# Glioma classification

classifyObj <-  reactive({
    inFile <- parseFilePaths(get.volumes(isolate({input$workingDir})), input$classifyObj)
    if (nrow(inFile) == 0) return(NULL)
    file <- as.character(inFile$datapath)
    withProgress(message = 'Reading file',
                 detail = 'This may take a while...', value = 0, {
                     if(grepl("\\.csv", file)){
                         ret <- read_csv(file)
                     } else {
                         ret <- get(load(file))
                     }
                 })
    return(ret)
})

observe({
    shinyFileChoose(input, 'classifyObj', roots=get.volumes(input$workingDir), session=session,
                    restrictions=system.file(package='base'), filetypes=c('', "rda","csv"))
})

observeEvent(input$idatClassify, {
    output$idattbl <- DT::renderDataTable({
        met <- classifyObj()
        if(is.null(met)) return(NULL)
        if(class(met)[1] == "RangedSummarizedExperiment") {
            met <- assay(met) %>% as.matrix  %>% t
        } else {
            met <- met %>% as.matrix  %>% t
        }

        df.all <- NULL
        models <- c("glioma.idhmut.model","glioma.gcimp.model")
        for(i in models){
            model <- get(i)
            # If it is a Summarized Experiment object

            # keep only probes used in the model
            aux <- met[,colnames(met) %in% model$coefnames]

            # This should not happen!
            if(any(apply(aux,2,function(x) all(is.na(x))))) {
                print("NA columns")
                aux[,apply(aux,2,function(x) all(is.na(x)))] <- 0.5 # runif(1, 0, 1) # pick a random number between 0 and 1
            }
            if(any(apply(aux,2,function(x) any(is.na(x))))) {
                print("NA values")
                colMedians <- colMedians(aux,na.rm = T)
                x <- which(is.na(aux),arr.ind = T)
                for(l in 1:nrow(x)){
                    aux[x[l,1],x[l,2]] <- colMedians[x[l,2]]
                }
            }


            pred <- predict(model, aux)
            df <- data.frame(samples = rownames(aux), groups.classified = pred)

            if(is.null(df.all)) {
                df.all <- df
            } else {
                df.all <- merge(df, df.all, by = "samples")
            }
        }
        colnames(df.all) <- c("samples",models)
        return(createTable(df.all))
    })
})

