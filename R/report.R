# Create Report
#' @import ReporteRs ggplot2
#' @importFrom UpSetR upset
create.report <- function(query, path = "report", system) {

    # pages: main.html encode.html tcga.html table.html figures.html
    pb <- txtProgressBar(min = 0, max = 5, style = 3)

    # creating footer
    footer <- pot( 'Code licensed under ', format = textProperties(color='gray') ) +
        pot('GPL-3', format = textProperties(color='#428bca'),
            hyperlink = 'https://gnu.org/licenses/gpl.html' ) +
        pot('.', format = textProperties(color='gray') )

    # creating menu
    menu <- BootstrapMenu (title = 'biOMICs report',link = "main.html")
    menu <- addLinkItem (menu, label = 'Graphs', 'figures.html')
    menu <- addLinkItem (menu, label = 'TCGA data summary', 'tcga.html')
    menu <- addLinkItem (menu, label = 'Encode data summary', 'encode.html')
    menu <- addLinkItem (menu, label = 'Results table', 'table.html')

    # MAIN --------------------------------------------------------------
    doc = bsdoc(title = 'biomics' )
    doc = addBootstrapMenu( doc, menu )
    doc = addFooter( doc, value = footer, par.properties = parCenter( padding = 2 ))
    doc <- main.report(doc,query,system)
    writeDoc(doc, file =  file.path(path,"main.html" ))
    setTxtProgressBar(pb, 1)


    #---------------- ENCODE file summary ----------------
    doc = bsdoc( title = 'biomics' )
    doc = addBootstrapMenu(doc,menu)
    doc <- encode.report(doc,query[query$database=="encode",])
    doc = addFooter(doc, value = footer, par.properties = parCenter( padding = 2 ))
    writeDoc( doc, file = file.path(path,"encode.html" ))
    setTxtProgressBar(pb, 3)

    # TCGA --------------------------------------------------------------
    doc = bsdoc( title = 'biomics' )
    doc = addBootstrapMenu(doc,menu)
    doc <- tcga.report(doc,query[query$database=="tcga",])
    doc = addFooter(doc, value = footer, par.properties = parCenter( padding = 2 ))
    writeDoc( doc, file = file.path(path,"tcga.html" ))
    setTxtProgressBar(pb, 2)

    # Graphs ---------------------------------------------------
    doc = bsdoc( title = 'biomics' )
    doc = addBootstrapMenu(doc,menu)
    doc <- figures.report(doc,query)
    doc = addFooter(doc, value = footer, par.properties = parCenter( padding = 2 ))
    writeDoc( doc, file = file.path(path,"figures.html" ))
    setTxtProgressBar(pb, 4)

    # Table -----------------------------------------------
    doc = bsdoc( title = 'biomics' )
    doc = addBootstrapMenu(doc, menu )
    doc = addFooter(doc, value = footer, par.properties = parCenter( padding = 2 ))
    doc = addFlexTable( doc, vanilla.table(query) )
    writeDoc(doc, file =  file.path(path,"table.html"))
    setTxtProgressBar(pb, 5)

    message(paste0("Report saved in folder:",path))
}

main.report <- function(doc,query,system){
    message("Creating main page page")
    mkd = "## Summary of results\n"
    mkd = paste0(mkd,paste0("Mapped to: **", system ,"**"))
    doc = addMarkdown( doc, text = mkd,
                       default.par.properties = parProperties(text.align = "justify",
                                                              padding.left = 0) )

    mkd = "### Number of terms in the system"
    doc = addMarkdown( doc, text = mkd,
                       default.par.properties = parProperties(text.align = "justify",
                                                              padding.left = 0) )
    doc = addPlot( doc, function() barplot(c(length(unique(query[query$database == "tcga","Sample"])),
                                             length(unique(query[query$database == "encode","Sample"])),
                                             length(unique(query[query$database == "roadmap","Sample"]))),
                                           xlab="Nuber of terms",
                                           legend.text = c(paste0("tcga (n=",length(unique(query[query$database == "tcga","Sample"])),")"),
                                                           paste0("encode (n=",length(unique(query[query$database == "encode","Sample"])),")"),
                                                           paste0("roadmap (n=",length(unique(query[query$database == "roadmap","Sample"])),")")),
                                           col = rainbow(3)),
                   width = 5, height = 5 )

    mkd = "### Number of samples in the system"

    doc = addMarkdown( doc, text = mkd,
                       default.par.properties = parProperties(text.align = "justify",
                                                              padding.left = 0) )
    doc = addPlot( doc, function() barplot(table(query$database), 	xlab="Nuber of samples",
                                           legend.text = c(paste0("tcga (n=",nrow(query[query$database == "tcga",]),")"),
                                                           paste0("encode (n=",nrow(query[query$database == "encode",]),")"),
                                                           paste0("roadmap (n=",nrow(query[query$database == "roadmap",]),")")),
                                           col = rainbow(3)),
                   width = 5, height = 5 )


    mkd = "## Downloading the samples"
    doc = addMarkdown( doc, text = mkd,
                       default.par.properties = parProperties(text.align = "justify",
                                                              padding.left = 0) )
    mkd = "The following subsections shows how the data seached can be downloaded using different biocondcutor packages"
    doc = addMarkdown( doc, text = mkd,
                       default.par.properties = parProperties(text.align = "justify",
                                                              padding.left = 0) )


    mkd = "### Downloading roadmap data with AnnotationHub\n"
    doc = addMarkdown( doc, text = mkd,
                       default.par.properties = parProperties(text.align = "justify",
                                                              padding.left = 0) )

    myrcode = 'library(pbapply)
library(AnnotationHub)
ah = AnnotationHub()
epiFiles <- sapply(unique(query[query$database=="roadmap",]$ID), function(x){names(query(ah,c(x,"EpigenomeRoadMap","narrowPeak","H3K4me2"))) })
epi.marks <- pblapply(unlist(epiFiles), function(x){ah[[x]]})'

    doc = addRScript(doc, text = myrcode )

    mkd = "### Downloading TCGA data with TCGAbiolinks\n"
    doc = addMarkdown( doc, text = mkd,
                       default.par.properties = parProperties(text.align = "justify",
                                                              padding.left = 0) )

    myrcode = 'library(tcgabiolinks)
tcga <- TCGAquery(tumor=query[query$database =="tcga",]$Sample,platform = query[query$database =="tcga",]$Experiment,level = 3)
TCGAdownload(tcga)'

    doc = addRScript(doc, text = myrcode )

    mkd = "### Downloading ENCODE data with biOMICs\n"
    doc = addMarkdown( doc, text = mkd,
                       default.par.properties = parProperties(text.align = "justify",
                                                              padding.left = 0) )

    myrcode = 'library(biOMICs)
encode <- query[query$database =="encode",]
biOmicsDownload(encode[1,],enc.file.type="bigWig")'

    doc = addRScript(doc, text = myrcode )
    message("Completed main page page")
    return(doc)
}

figures.report <-  function(doc, query){
    message("Creating figures summary page")
    report.plot <- function(data, col, title){
        .e <- environment()

        g <- ggplot(data, aes(factor(database),
                              fill = data[,col]),
                    environment = .e) +
            geom_bar(position = "fill") +
            theme_bw() +
            theme(panel.border = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.line = element_line(colour = "black"),
                  legend.key = element_rect(colour = 'white'),
                  legend.justification=c(1,1), legend.position=c(1,1)) +
            ggtitle(title) +
            labs(x="Database", y="Percentage of samples") +
            scale_fill_discrete(name=col) +
            theme(legend.direction ="vertical",legend.position = "bottom")+
            guides(fill=guide_legend(ncol=3))

        return(g)
    }

    mkd = "# Summary plots"
    doc = addMarkdown(doc, text = mkd, default.par.properties = parProperties(text.align = "justify", padding.left = 0) )

    # add a plot into doc
    doc = addPlot(doc, function() plot(report.plot(query,"Experiment","Experiment per database")), width = 15, height = 10 )

    doc = addPlot(doc, function() plot(report.plot(query,"Sample","Samples per database")), width = 15, height = 10 )

    mkd = "# Summary plots by platform "
    doc = addMarkdown( doc, text = mkd,
                       default.par.properties = parProperties(text.align = "justify", padding.left = 0) )

    for(i in unique(platforms$Standard)){
        aux <- subset(query,subset = query$Experiment %in%
                          platforms[platforms$Standard == i,"Platform"])
        if(nrow(aux) >0) {
            mkd = paste0("## Summary for ",i)
            doc = addMarkdown( doc, text = mkd,
                               default.par.properties = parProperties(text.align = "justify",
                                                                      padding.left = 0) )


            # add a plot into doc
            doc = addPlot( doc, function() plot(report.plot(aux,"Experiment",paste0("Experiment per database - ",i))),
                           width = 4 *  length(unique(aux$database)), height = 5 )

            # add a plot into doc
            doc = addPlot( doc, function() plot(report.plot(aux,"Sample",paste0("Samples per database - ",i))),
                           width = 4 *  length(unique(aux$database)), height = 5 + length(unique(aux$Sample))/10 )

        }
    }
    message("Completed figures summary page")
    return(doc)
}
#' @importFrom UpSetR upset
tcga.report <- function(doc, query){
    message("Creating TCGA summary page")
    mkd = "# Summary plots"
    doc = addMarkdown(doc, text = mkd, default.par.properties = parProperties(text.align = "justify", padding.left = 0) )

    for(level in 1:3){
        mkd = paste0("## Summary for level ",level, " data")
        doc = addMarkdown( doc, text = mkd,
                           default.par.properties = parProperties(text.align = "justify", padding.left = 0) )
        x <- TCGAquery(unique(query$Sample),level=level)
        patient <- unique(substr(unlist(stringr::str_split(x$barcode,",")),1,15))
        platform <- unique(x$Platform)
        df <- as.data.frame(matrix(0,nrow=length(patient),ncol=length(platform)+1))
        colnames(df) <- c("patient", platform)
        df$patient <- patient
        for (i in patient){
            idx <- grep(i,x$barcode)
            plat <- x[idx,"Platform"]
            for (j in plat){
                df[df$patient == i,j] <- 1
            }
        }

        doc = addPlot(doc,
                      function() {
                          table.code <- c("Primary solid Tumor","Recurrent Solid Tumor",
                                          "Primary Blood Derived Cancer - Peripheral Blood",
                                          "Recurrent Blood Derived Cancer - Bone Marrow",
                                          "Additional - New Primary",
                                          "Metastatic","Additional Metastatic",
                                          "Human Tumor Original Cells",
                                          "Primary Blood Derived Cancer - Bone Marrow",
                                          "Blood Derived Normal","Solid Tissue Normal",
                                          "Buccal Cell Normal","EBV Immortalized Normal",
                                          "Bone Marrow Normal","Control Analyte",
                                          "Recurrent Blood Derived Cancer - Peripheral Blood",
                                          "Cell Lines","Primary Xenograft Tissue",
                                          "Cell Line Derived Xenograft Tissue")

                          names(table.code) <- c('01','02','03','04','05','06','07','08','09','10',
                                                 '11','12','13','14','20','40','50','60','61')
                          tab <- table(substr(df$patient,14,15))
                          names(tab) <- table.code[names(tab)]
                          tab <- tab[!is.na(names(tab))]
                          tab <- (as.data.frame(tab))
                          tab$type <- rownames(tab)
                          colnames(tab) <- c("Freq","type")
                          p <- ggplot(tab, aes(x=type,y = Freq,fill=type)) + geom_bar(stat="identity") +
                              theme_bw() +theme(panel.border = element_blank(),
                                                panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(),
                                                axis.line = element_line(colour = "black"),
                                                legend.key = element_rect(colour = 'white'),
                                                legend.justification=c(1,1),
                                                legend.position=c(1,1),
                                                text = element_text(size=16),
                                                axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Type of sample") +
                              scale_fill_brewer(palette="Set1") + guides(fill=FALSE) + geom_text(aes(label = Freq), size = 4)
                          plot(p)

                      },width = 6, height = 5 )

        df$patient <- NULL

        doc = addPlot(doc,
                      function() {
                          upset(df, nsets = length(platform),
                                number.angles = 30,
                                nintersects = 100,
                                point.size = 3, name.size = 12,
                                line.size = 1,
                                mainbar.y.label = "Platform Intersections",
                                sets.x.label = "Samples Per Platform",
                                order.by = "freq",
                                decreasing = T,
                                #group.by = "sets",
                                sets.bar.color = "#56B4E9")

                      },width = 20, height = 15 )
    }
    message("Completed TCGA summary page")
    return(doc)
}


#' @importFrom UpSetR upset
encode.report <- function(doc, query){
    message("Creating encode summary page")
    mkd = paste0("## Summary for encode")
    doc = addMarkdown( doc, text = mkd,
                       default.par.properties = parProperties(text.align = "justify", padding.left = 0) )


    encode.db.files$dataset <- gsub("experiments","",encode.db.files$dataset)
    encode.db.files$dataset <- gsub("/","",encode.db.files$dataset)
    encode.db.files <- encode.db.files[encode.db.files$dataset %in% query$ID,]
    encode.tab <- merge(encode.db,encode.db.files[,1:4],by.x="accession",by.y="dataset")

    for(j in unique(encode.db$organism)) {
        mkd = paste0("### Organism:", j)
        doc = addMarkdown(doc, text = mkd,
                          default.par.properties = parProperties(text.align = "justify", padding.left = 0) )

        tab <- encode.tab[encode.tab$organism == j,]
        tab$files <- paste(tab$output_type, tab$title, tab$file_format, sep=",")

        for(i in unique(tab$accession)){

            tab[tab$accession == i,"files"] <- paste(tab[tab$accession == i, "files"], collapse = "\n")
        }
        tab <- tab[!duplicated(tab$accession),]

        tab <- tab[,-c(8:10)]
        colnames(tab)[8] <- "Files: description,title,type"

        options( "ReporteRs-fontsize" = 10 )
        tab <- vanilla.table(tab)


        tab = setFlexTableWidths( tab, widths = c(1, # accession
                                                  1, # biosample
                                                  2, # assay
                                                  1, # lab
                                                  1, # target
                                                  2, # description
                                                  1,   # organism
                                                  5))  # files
        tab = setZebraStyle( tab, odd = "#E1EEf4", even = "white" )
        doc = addFlexTable( doc, tab,width=50)
    }
    message("Completed encode summary page")

    return(doc)
}



