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

observeEvent(input$gliomaClassify, {
    closeAlert(session,"gliomaAlert")
    output$gliomatbl <- DT::renderDataTable({
        met <- classifyObj()
        if(is.null(met)) {
            createAlert(session, "gliomaAlert", "gliomamessage", title = "Error", style =  "danger",
                        content = "Please select a file", append = FALSE)
            return(NULL)
        }
        if(class(met)[1] == "RangedSummarizedExperiment") {
            met <- assay(met) %>% as.matrix  %>% t
        } else {
            met <- met %>% as.matrix  %>% t
        }

        withProgress(message = 'Classifying samples',
                     detail = 'Selecting models', value = 0, {
                         tryCatch({
                             df.all <- NULL
                             models <- c("idh","gcimp","idhwt","idhmut")
                             models <- paste("glioma",models,"model",sep = ".")
                             data(list = models)
                             for(i in models){
                                 incProgress(0.2, message = "Classifying into glioma subtypes",detail = paste0("Model: ", models[i]))
                                 model <- get(i)
                                 # If it is a Summarized Experiment object

                                 # keep only probes used in the model
                                 aux <- met[,colnames(met) %in% colnames(model$trainingData),drop=FALSE]

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
                                 df <- data.frame(samples = rownames(aux), groups.classified = pred,stringsAsFactors = FALSE)

                                 if(is.null(df.all)) {
                                     df.all <- df
                                 } else {
                                     df.all <- merge(df.all,df, by = "samples")
                                 }
                             }
                             incProgress(0.2, message = "Classifying into glioma subtypes",detail = paste0("Preparing final table"))
                             colnames(df.all) <- c("samples",models)
                             fctr.cols <- sapply(df.all, is.factor)
                             df.all[, fctr.cols] <- sapply(df.all[, fctr.cols], as.character)
                             df.all[grep("6|5|4",df.all$glioma.idh.model),c("glioma.gcimp.model","glioma.idhmut.model")]  <- NA
                             df.all[grep("3|2|1",df.all$glioma.idh.model),c("glioma.idhwt.model")]  <- NA
                             df.all[grep("3",df.all$glioma.idhmut.model),c("glioma.gcimp.model")]  <- "Codel"
                             df.all[grep("1",df.all$glioma.idhwt.model),c("glioma.idhwt.model")]  <- "Classic-like"
                             df.all[grep("2",df.all$glioma.idhwt.model),c("glioma.idhwt.model")]  <- "Mesenchymal-like"
                             df.all[grep("3",df.all$glioma.idhwt.model),c("glioma.idhwt.model")]  <- "PA-like"
                             # Final column with results
                             df.all$glioma.DNAmethylation.subtype <- NA
                             df.all$glioma.DNAmethylation.subtype <- df.all$glioma.idhwt.model
                             idx <- which(is.na(df.all$glioma.DNAmethylation.subtype))
                             df.all$glioma.DNAmethylation.subtype[idx] <- df.all$glioma.gcimp.model[idx]
                         }, error = function(e){
                             createAlert(session, "gliomaAlert", "gliomamessage", title = "Error", style =  "error",
                                         content =  paste0("Erro: ", "<br><ul>", paste(e, collapse = "</ul><ul>"),"</ul>"), append = FALSE)
                         })
                     })
        return(createTable(df.all,"Glioma classification"))
    })
})

