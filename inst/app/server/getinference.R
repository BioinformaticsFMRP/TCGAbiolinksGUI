# -------------------------------------------------
# Network inference
# -------------------------------------------------
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# File selection
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observe({
    shinyFileChoose(input, 'networkfile', roots=get.volumes(input$workingDir),
                    session=session, restrictions=system.file(package='base'),
                    filetypes=c('rda', 'RData'))
})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
# Network inference
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
observeEvent(input$networkBt , {
    print("Start inference")
        inFile <- input$networkfile

        load(paste0("~",paste(unlist(inFile$files), collapse="/")))
 
        data <- t(assays(data)[[2]])

        myadj <- TCGAanalyze_networkInference(data, optionMethod=input$networkInferenceMethodRb)
        
        out.filename <- paste0(paste("network_results_method", input$networkInferenceMethodRb,sep="_"),".rda")
                print(out.filename)

        getPath <- parseDirPath(get.volumes(isolate({input$workingDir})), input$workingDir)
        if (length(getPath) == 0) getPath <- paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI")

        save(myadj, file = out.filename)
        print("file written")

})


