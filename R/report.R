# Create Report
#' @import ReporteRs ggplot2
#' @importFrom UpSetR upset
create.report <- function(query, path = "report", system) {

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

    # creating footer
    footer <- pot( 'Code licensed under ', format = textProperties(color='gray') ) +
        pot('GPL-3', format = textProperties(color='#428bca'),
            hyperlink = 'https://gnu.org/licenses/gpl.html' ) +
        pot('.', format = textProperties(color='gray') )

    # creating menu
    menu <- BootstrapMenu (title = 'biOMICs report',link = "main.html")
    menu <- addLinkItem (menu, label = 'Graphs', 'figures.html')
    menu <- addLinkItem (menu, label = 'TCGA data summary', 'tcga.html')
    menu <- addLinkItem (menu, label = 'Table', 'table.html')

    message("Creating report....")

    # MAIN --------------------------------------------------------------
    doc = bsdoc(title = 'biomics' )
    doc = addBootstrapMenu( doc, menu )
    doc = addFooter( doc, value = footer, par.properties = parCenter( padding = 2 ))

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


    mkd = "## Downloading the samples
The next subsections shows how the data seached can be downloaded using different biocondcutor packages"
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


    # write the doc
    writeDoc( doc, file =  file.path(path,"main.html" ))

    # TGCA data summary
    # Creation of doc, a docx object
    doc = bsdoc( title = 'biomics' )
    doc = addBootstrapMenu(doc,menu)

    mkd = "# Summary plots"
    doc = addMarkdown(doc, text = mkd, default.par.properties = parProperties(text.align = "justify", padding.left = 0) )
    doc <- tcga.report(doc,query[query$database=="tcga",])

    doc = addFooter(doc, value = footer, par.properties = parCenter( padding = 2 ))

    # write the doc
    writeDoc( doc, file = file.path(path,"tcga.html" ))

    # Graphs ---------------------------------------------------

    # Creation of doc, a docx object
    doc = bsdoc( title = 'biomics' )
    doc = addBootstrapMenu(doc,menu)

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

    doc = addFooter(doc, value = footer, par.properties = parCenter( padding = 2 ))

    # write the doc
    writeDoc( doc, file = file.path(path,"figures.html" ))

    # Table -----------------------------------------------
    doc = bsdoc( title = 'biomics' )
    doc = addBootstrapMenu(doc, menu )
    doc = addFooter(doc, value = footer, par.properties = parCenter( padding = 2 ))

    # add into doc first 10 lines of iris
    doc = addFlexTable( doc, vanilla.table(query) )
    # write the doc
    writeDoc(doc, file =  file.path(path,"table.html"))

    message(paste0("Report saved in folder:",path))
}

#' @importFrom UpSetR upset
tcga.report <- function(doc, query){

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
    return(doc)
}


