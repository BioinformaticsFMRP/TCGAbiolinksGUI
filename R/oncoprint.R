#' Creating a oncoprint
#' @param mut Mutation file (see TCGAquery_maf from TCGAbiolinks)
#' @param genes Gene list
#' @param filename name of the pdf
#' @param color named vector for the plot
#' @param height pdf height
#' @param rm.empty.columns If there is no alteration in that sample, whether remove it on the oncoprint
#' @param show.row.barplot  Show barplot annotation on rows?
#' @param show.column.names Show column names? Default: FALSE
#' @param font.size Size of the fonts
#' @param annotation Matrix or data frame with the annotation.
#' Should have a column bcr_patient_barcode with the same ID of the mutation object
#' @param annotation.position Position of the annotation "bottom" or "top"
#' @param label.title Title of the label
#' @importFrom ComplexHeatmap oncoPrint draw HeatmapAnnotation
#' @importFrom grid gpar grid.rect
#' @importFrom reshape2 dcast acast
#' @examples
#' mut <- TCGAbiolinks::TCGAquery_maf(tumor = "GBM",
#'                                    archive.name = "ucsc.edu_GBM.IlluminaGA_DNASeq_automated.Level_2.1.1.0")
#' create.oncoprint(mut = mut, genes = mut$Hugo_Symbol[1:10])
#' create.oncoprint(mut = mut, genes = mut$Hugo_Symbol[1:10],
#'                  filename = "onco.pdf",
#'                  color=c("DEL"="purple","INS"="yellow","SNP"="brown"))
#' clin <- TCGAbiolinks::TCGAquery_clinic("gbm","clinical_patient")
#' clin <- clin[,c("bcr_patient_barcode","disease","gender")]
#' create.oncoprint(mut = mut, genes = mut$Hugo_Symbol[1:10],
#'                  filename = "onco.pdf",
#'                  annotation = clin,
#'                  color=c("DEL"="purple","INS"="yellow","SNP"="brown"),font.size=10)
#'
#' @export
create.oncoprint <- function (mut,
                              genes,
                              filename,
                              color,
                              annotation.position = "bottom",
                              annotation,
                              height,
                              rm.empty.columns = FALSE,
                              show.column.names = FALSE,
                              show.row.barplot = TRUE,
                              label.title = "Mutation",
                              font.size = 16){


    if(missing(mut))   stop("Missing mut argument")
    mut$value <- 1
    if(!rm.empty.columns){
        mat <- dcast(mut, Tumor_Sample_Barcode + Hugo_Symbol ~ Variant_Type,value.var = "value",fill = 0,drop = FALSE)
    } else {
        mat <- dcast(mut, Tumor_Sample_Barcode + Hugo_Symbol ~ Variant_Type,value.var = "value",fill = 0,drop = TRUE)
    }
    if(!missing(genes) & !is.null(genes)) mat <- subset(mat, mat$Hugo_Symbol %in% genes)

    # mutation in the file
    columns <- colnames(mat)[-c(1:2)]

    # value will be a collum with all the mutations
    mat$value <- ""
    for ( i in columns){
        mat[,i] <-  replace(mat[,i],mat[,i]>0,paste0(i,";"))
        mat[,i] <-  replace(mat[,i],mat[,i]==0,"")
        mat$value <- paste0(mat$value,mat[,i])
    }

    # now we have a matrix with pairs samples/genes mutations
    # we want a matrix with samples vs genes mutations with the content being the value
    mat <- acast(mat, Tumor_Sample_Barcode~Hugo_Symbol, value.var="value",fill="")
    rownames(mat) <-  substr(rownames(mat),1,12)

    alter_fun = list(
        background = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
        },
        INS = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = color["INS"], col = NA))
        },
        DEL = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = color["DEL"], col = NA))
        },
        SNP = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = color["SNP"], col = NA))
        },
        DNP = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = color["DNP"], col = NA))
        }
    )

    # After the gene selection, some of the mutation might not exist
    # we will remove them to make the oncoprint work
    mutation.type <- c()
    for (i in columns){
        if(length(grep(i,mat)) > 0) mutation.type <- c(mutation.type,i)
    }
    # get only the colors to the mutations
    # otherwise it gives errors

    if(missing(color)){
        color <- c(rainbow(length(mutation.type)))
        names(color) <- mutation.type
    } else{
        color <- color[mutation.type]
    }
    alt = intersect(names(alter_fun), c("background",as.character(mutation.type)))
    alter_fun <- alter_fun[alt]

    # header are samples, rows genes
    mat <- t(mat)

    if(!missing(height)) height <- length(genes)/2
    if(!missing(filename)) pdf(filename,width = 20,height = height)

    if(missing(annotation)) annotation <- NULL
    if(!is.null(annotation)){
        idx <- match(substr(colnames(mat),1,12),annotation$bcr_patient_barcode)

        stopifnot(all(annotation$bcr_patient_barcode[idx] ==substr(colnames(mat),1,12) ))
        annotation <- annotation[idx,]

        annotation$bcr_patient_barcode <- NULL

        n.col <- sum(sapply(colnames(annotation), function(x) {
            length(unique(annotation[,x]))
        }))

        # add automatic colors: not working

        get.color <- function(df,col){
            idx <- which(colnames(df) == col)
            start <- 1
            if(idx != 1) start <- length(unique(unlist(c(df[,1:(idx-1)])))) + 1
            end <- start + length(unique(df[,col])) -1
            diff.colors <- c("dimgray","thistle","deeppink3","magenta4","lightsteelblue1","black","chartreuse","lightgreen","maroon4","darkslategray","lightyellow3","darkslateblue","firebrick1","aquamarine","dodgerblue4","bisque4","moccasin","indianred1","yellow","gray93","cyan","darkseagreen4","lightgoldenrodyellow","lightpink","sienna1","darkred","palevioletred","tomato4","blue","mediumorchid4","royalblue1","magenta2","darkgoldenrod1")
            return(diff.colors[start:end])
        }
        col.annot <- lapply(colnames(annotation), function(x) {
            #idx <- which(colnames(annotation) == x) - 1
            #print(idx/n.col)
            ret <- get.color(annotation,x)
            #ret <- rainbow(length(unique(annotation[,x])),start = idx/n.col,alpha=0.5)
            names(ret) <- as.character(unique(annotation[,x]))
            return(ret)
        })
        names(col.annot) <-  colnames(annotation)

        annotHeatmap <- HeatmapAnnotation(df=annotation,
                                          col=col.annot,
                                          annotation_legend_param=list(title_gp=gpar(fontsize=font.size,
                                                                                     fontface="bold"),
                                                                       labels_gp=gpar(fontsize=font.size),#sizelabels
                                                                       grid_height=unit(8,"mm")))
    }
    if(is.null(annotation)){
        p <- oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                       row_order = NULL,
                       remove_empty_columns = FALSE,
                       show_column_names = show.column.names,
                       show_row_barplot = show.row.barplot,
                       column_order = NULL, # Do not sort the columns
                       alter_fun = alter_fun, col = color,
                       row_names_gp = gpar(fontsize = font.size),  # set size for row names
                       pct_gp = gpar(fontsize = font.size), # set size for percentage labels
                       axis_gp = gpar(fontsize = font.size),# size of axis
                       #column_title = "OncoPrint for TCGA LGG, genes in Glioma signaling",
                       #column_title_gp = gpar(fontsize = 11),
                       row_barplot_width = unit(2, "cm"), #size barplot
                       heatmap_legend_param = list(title = label.title, at = names(color),
                                                   labels = names(color),
                                                   title_gp = gpar(fontsize = font.size, fontface = "bold"),
                                                   labels_gp = gpar(fontsize = font.size), # size labels
                                                   grid_height = unit(8, "mm") # vertical distance labels
                       )
        )
    } else if(!is.null(annotation) & annotation.position == "bottom"){

        p <- oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                       row_order = NULL,
                       remove_empty_columns = FALSE,
                       show_row_barplot = show.row.barplot,
                       show_column_names = show.column.names,
                       column_order = NULL, # Do not sort the columns
                       alter_fun = alter_fun, col = color,
                       row_names_gp = gpar(fontsize = font.size),  # set size for row names
                       pct_gp = gpar(fontsize = font.size), # set size for percentage labels
                       axis_gp = gpar(fontsize = font.size),# size of axis
                       #column_title = "OncoPrint for TCGA LGG, genes in Glioma signaling",
                       #column_title_gp = gpar(fontsize = 11),
                       row_barplot_width = unit(2, "cm"), #size barplot
                       bottom_annotation = annotHeatmap,
                       heatmap_legend_param = list(title = label.title, at = names(color),
                                                   labels = names(color),
                                                   title_gp = gpar(fontsize = font.size, fontface = "bold"),
                                                   labels_gp = gpar(fontsize = font.size), # size labels
                                                   grid_height = unit(8, "mm") # vertical distance labels
                       )
        )

    } else {
        p <- oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                       row_order = NULL,
                       remove_empty_columns = FALSE,
                       show_column_names = show.column.names,
                       show_row_barplot = show.row.barplot,
                       column_order = NULL, # Do not sort the columns
                       alter_fun = alter_fun, col = color,
                       row_names_gp = gpar(fontsize = font.size),  # set size for row names
                       pct_gp = gpar(fontsize = font.size), # set size for percentage labels
                       axis_gp = gpar(fontsize = font.size),# size of axis
                       #column_title = "OncoPrint for TCGA LGG, genes in Glioma signaling",
                       #column_title_gp = gpar(fontsize = 11),
                       row_barplot_width = unit(2, "cm"), #size barplot
                       top_annotation = annotHeatmap,
                       heatmap_legend_param = list(title = label.title, at = names(color),
                                                   labels = names(color),
                                                   title_gp = gpar(fontsize = font.size, fontface = "bold"),
                                                   labels_gp = gpar(fontsize = font.size), # size labels
                                                   grid_height = unit(8, "mm") # vertical distance labels
                       )
        )
    }

    if(!missing(filename)) dev.off()
    draw(p)
}
