# Subtype
observeEvent(input$tcgaSubtypeBt, {
    updateCollapse(session, "collapseTCGASubtype", open = "Subtype data: results")
    output$tcgaSubtypetbl <- renderDataTable({
        tbl <- subtype.result()
        tumor <- isolate({input$tcgasubtypeFilter})
        if(!is.null(tbl)) {
            closeAlert(session, "tcgaSubtypeAlert")
            doi <- c("acc"="Comprehensive Pan-Genomic Characterization of Adrenocortical Carcinoma<br>doi:10.1016/j.ccell.2016.04.002",
                     "aml"="Genomic and epigenomic landscapes of adult de novo acute<br>doi:10.1056/NEJMoa1301689",
                     "blca"="Comprehensive molecular characterization of urothelial bladder <br>doi:10.1038/nature12965",
                     "brca"="Comprehensive molecular portraits of human breast tumours<br>doi:10.1038/nature11412",
                     "coad"="Comprehensive molecular characterization of human colon<br>doi:10.1038/nature11252",
                     "gbm"="Molecular Profiling Reveals Biologically Discrete Subsets and Pathways of Progression in Diffuse Glioma.<br>doi:10.1016/j.cell.2015.12.028",
                     "lgg"="Molecular Profiling Reveals Biologically Discrete Subsets and Pathways of Progression in Diffuse Glioma.<br>doi:10.1016/j.cell.2015.12.028",
                     "hnsc"="Comprehensive genomic characterization of head and neck<br>doi:10.1038/nature14129",
                     "kich"="The somatic genomic landscape of chromophobe renal cell carcinoma<br>doi:10.1016/j.ccr.2014.07.014",
                     "kirc"="Comprehensive molecular characterization of clear cell renal cell carcinoma<br>doi:10.1038/nature12222",
                     "kirp"="Comprehensive Molecular Characterization of Papillary Renal-Cell Carcinoma<br>doi:10.1056/NEJMoa1505917",
                     "lihc"="",
                     "luad"="Comprehensive molecular profiling of lung adenocarcinoma <br>doi:10.1038/nature13385",
                     "lusc"="Comprehensive genomic characterization of squamous cell lung cancers<br>doi:10.1038/nature11404",
                     "ovca"= "Integrated genomic analyses of ovarian carcinoma<br>doi:10.1038/nature10166",
                     "pancan"="Multiplatform analysis of 12 cancer types reveals molecular <br>doi:10.1016/j.cell.2014.06.049",
                     "pcpg"="Comprehensive Molecular Characterization of Pheochromocytoma and Paraganglioma<br>http://dx.doi.org/10.1016/j.ccell.2017.01.001",
                     "prad"="The Molecular Taxonomy of Primary Prostate Cancer<br>doi:10.1016/j.cell.2015.10.025",
                     "skcm"="Genomic Classification of Cutaneous Melanoma<br>doi:10.1016/j.cell.2015.05.044",
                     "stad"="Comprehensive molecular characterization of gastric adenocarcinoma<br>doi:10.1038/nature13480",
                     "thca"="Integrated genomic characterization of papillary thyroid carcinoma<br>doi:10.1016/j.cell.2014.09.050",
                     "ucec"="Integrated genomic characterization of endometrial carcinoma<br>doi:10.1038/nature12113",
                     "ucs"="")
            if (isolate({input$saveSubtypeRda}) || isolate({input$saveSubtypeCsv})) {
                save.message <- ""
                getPath <- parseDirPath(get.volumes(isolate({input$workingDir})), input$workingDir)
                if (length(getPath) == 0) getPath <- paste0(Sys.getenv("HOME"),"/TCGAbiolinksGUI")
                filename <- file.path(getPath,paste0(tumor,"_subtype.rda"))
                if (isolate({input$saveSubtypeRda})) {
                    save(tbl, file = filename)
                    save.message <- paste0(save.message,"<br> File created: ", filename)
                }
                if (isolate({input$saveSubtypeCsv})) {
                    write_csv(tbl, path = gsub("rda","csv",filename))
                    save.message <- paste0(save.message,"<br> File created: ", gsub("rda","csv",filename))
                }
                createAlert(session, "tcgaSubtypemessage", "tcgaSubtypeAlert", title = paste0("Success"), style =  "success",
                            content = paste0(save.message,"<br>Source of the data: ", doi[tumor]), append = TRUE)

            } else {
                createAlert(session, "tcgaSubtypemessage", "tcgaSubtypeAlert", title = "Source of the data", style =  "success",
                            content = paste0(doi[tumor]), append = TRUE)
            }

            return(tbl)
        }

    },
    options = list(pageLength = 10,
                   scrollX = TRUE,
                   jQueryUI = TRUE,
                   pagingType = "full",
                   lengthMenu = list(c(10, 20, -1), c('10', '20', 'All')),
                   language.emptyTable = "No results found",
                   "dom" = 'T<"clear">lfrtip',
                   "oTableTools" = list(
                       "sSelectedClass" = "selected",
                       "sRowSelect" = "os",
                       "sSwfPath" = paste0("//cdn.datatables.net/tabletools/2.2.4/swf/copy_csv_xls.swf"),
                       "aButtons" = list(
                           list("sExtends" = "collection",
                                "sButtonText" = "Save",
                                "aButtons" = c("csv","xls")
                           )
                       )
                   )
    ), callback = "function(table) {table.on('click.dt', 'tr', function() {Shiny.onInputChange('allRows',table.rows('.selected').data().toArray());});}"
    )
})

subtype.result <-  reactive({
    tumor <- isolate({input$tcgasubtypeFilter})
    tbl <- data.frame()

    result = tryCatch({
        tbl <- rbind(tbl, TCGAquery_subtype(tumor = tumor))

    }, error = function(e) {
        createAlert(session, "tcgaSubtypemessage", "tcgaSubtypeAlert", title = "No results found", style =  "warning",
                    content = "Sorry there are subtypes for your query.", append = FALSE)
    })

    if(is.null(tbl)){
        createAlert(session, "tcgaSubtypemessage", "tcgaSubtypeAlert", title = "No results found", style =  "warning",
                    content = "Sorry there are subtypes for your query.", append = FALSE)
        return()
    } else if(nrow(tbl) ==0) {
        createAlert(session, "tcgaSubtypemessage", "tcgaSubtypeAlert", title = "No results found", style =  "warning",
                    content = "Sorry there are subtypes for your query.", append = FALSE)
        return()
    }
    return(tbl)
})

observeEvent(input$subtypePlotCol, {
    col <- isolate({input$subtypePlotCol})
    if(!is.null(col) & str_length(col) > 1) {
        updateCollapse(session, "collapseTCGASubtype", open = "Subtype data: Summary")
    }
    output$subtypeview <- renderPlotly({
        tbl <- subtype.result()
        if(is.null(tbl)) return(plotly_empty())
        if(is.null(col)) return(plotly_empty())

        if(str_length(col) < 2) return(plotly_empty())
        df <- as.data.frame(dplyr::count_(tbl, eval(col)))
        colnames(df) <- c("Var1","Freq")
        df$Var1 <- as.character(df$Var1)
        df$Var1[is.na(df$Var1)] <- "No value"
        save(df,file = "test.rda")
        p <- plot_ly(df, labels = ~Var1, values = ~Freq, type = 'pie',
                     textposition = 'inside',
                     textinfo = 'label+percent',
                     insidetextfont = list(color = '#FFFFFF'),
                     hoverinfo = 'text',
                     text=~paste0(Var1,"\n",Freq),
                     marker = list(colors = colors,
                                   line = list(color = '#FFFFFF', width = 1)),
                     showlegend = FALSE) %>%
            layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                   yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)) %>%
            config(displayModeBar = F)
    })
})
observe({
    updateSelectizeInput(session, 'subtypePlotCol', choices =  colnames(subtype.result()), server = TRUE)
})
