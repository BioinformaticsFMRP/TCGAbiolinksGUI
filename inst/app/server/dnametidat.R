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
    IDATs <- dir(path = baseDir,full.names = T,recursive = T,pattern = "idat")
    IDATprefixes <- unique(gsub("_Grn.idat|_Red.idat","",IDATs))
    withProgress(message = 'Processing IDAT files with sesame',
                 detail = 'This may take a while...', value = 0, {
                     idat <- tryCatch({
                         betas <- openSesame(IDATprefixes,
                                    quality.mask = FALSE,
                                    nondetection.mask = FALSE,
                                    mask.use.tcga = FALSE,
                                    pval.threshold = 0.05)
                         colnames(betas) <- basename(IDATprefixes)
                         betas <- TCGAbiolinks:::makeSEFromDNAMethylationMatrix(betas,genome = input$IDATgenome,met.platform = input$IDATplatform)
                         betas
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
        inFile <- input$IDATfolder
        if(class(inFile) != "list") return(NULL)
        baseDir  <- parseDirPath(get.volumes(input$workingDir), input$IDATfolder)
        IDATs <- dir(path = baseDir,full.names = T,recursive = T,pattern = "idat")
        IDATprefixes <- unique(gsub("_Grn.idat|_Red.idat","",IDATs))
        df <- data.frame(samples = IDATprefixes)
        createTable(df)
    })
})

observeEvent(input$IDATplatform, {
    updateTextInput(session, "idatfilename", label = "Filename", value = gsub("450K|EPIC|27K",input$IDATplatform,input$idatfilename))
})
observeEvent(input$IDATgenome, {
    updateTextInput(session, "idatfilename", label = "Filename", value = gsub("hg38|hg19",input$IDATgenome,input$idatfilename))
})

observeEvent(input$idatnormalize, {
    closeAlert(session,"idatmessage")
    betas <- idat()
    if (is.null(idat)) return(NULL)

    withProgress(message = 'Normalizing IDAT files',
                 detail = 'This may take a while...', value = 0, {

                     fname <- isolate({input$idatfilename})
                     getPath <- parseDirPath(get.volumes(isolate({input$workingDir})), isolate({input$workingDir}))
                     if (length(getPath) == 0) getPath <- paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI")
                     fname <- file.path(getPath,fname)

                     incProgress(1, message = "Saving",detail = paste0("As: ", fname))
                     if(tools::file_ext(fname) == "csv"){
                         write.csv(assay(betas),file = fname,row.names = T)
                     } else  if(tolower(tools::file_ext(fname))  %in% c("rda","rdata")){
                         save(betas,file = fname)
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
