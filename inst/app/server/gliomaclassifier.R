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

        df.all <- NULL
        models <- c("idh","gcimp","idhwt","idhmut")
        models <- paste("glioma",models,"model",sep = ".")
        for(i in models){
            model <- get(i)
            # If it is a Summarized Experiment object

            # keep only probes used in the model
            aux <- met[,colnames(met) %in% colnames(model$trainingData)]

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
                df.all <- merge(df.all,df, by = "samples")
            }
        }
        colnames(df.all) <- c("samples",models)
        df.all[grep("6|5|4",df.all$glioma.idh.model),c("glioma.gcimp.model","glioma.idhmut.model")]  <- NA
        df.all[grep("3|2|1",df.all$glioma.idh.model),c("glioma.idhwt.model")]  <- NA
        df.all[grep("3",df.all$glioma.idhmut.model),c("glioma.gcimp.model")]  <- "Codel"

        return(createTable(df.all))
    })
})

