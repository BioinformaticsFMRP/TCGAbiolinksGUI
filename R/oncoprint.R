#' Creating a oncoprint
#' @param mut Mutation file (see TCGAquery_maf from TCGAbiolinks)
#' @param genes Gene list
#' @param filename name of the pdf
#' @param color named vector for the plot
#' @param height pdf height
#' @importFrom ComplexHeatmap oncoPrint
#' @importFrom grid gpar
#' @importFrom reshape2 dcast acast
#' @examples
#' mut <- TCGAbiolinks::TCGAquery_maf(tumor = "GBM", archive.name = "ucsc.edu_GBM.IlluminaGA_DNASeq_automated.Level_2.1.1.0")
#' create.oncoprint(mut = mut, genes = mut$Hugo_Symbol[1:10])
#' create.oncoprint(mut = mut, genes = mut$Hugo_Symbol[1:10],
#'                  filename = "onco.pdf",
#'                  color=c("DEL"="purple","INS"="yellow","SNP"="brown"))
#' @export
create.oncoprint <- function (mut,
                              genes,
                              filename,
                              color,
                              bottom_annotation,
                              height){

    if(missing(genes)) stop("Missing genes argument")
    if(missing(mut))   stop("Missing mut argument")
    if(missing(color))   color = c("DEL" = "#008000", "INS" = "red", "SNP" = "blue")

    mut <- subset(mut,mut$Hugo_Symbol %in% genes)
    mat <- dcast(mut, Tumor_Sample_Barcode + Hugo_Symbol ~ Variant_Type)

    columns <- c("DEL","INS","SNP")
    columns <- columns[columns %in% colnames(mat)]
    for( i in columns){
        mat[,i] <-  replace(mat[,i],mat[,i]>0,paste0(i,";"))
        mat[,i] <-  replace(mat[,i],mat[,i]==0,"")
    }

    mat$value <- paste0(mat$SNP,mat$INS,mat$DEL)
    mat <- acast(mat, Tumor_Sample_Barcode~Hugo_Symbol, value.var="value",fill="")
    rownames(mat) <-  substr(rownames(mat),1,12)

    alter_fun_list = list(
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
        }
    )

    mat <- t(mat)

    if(!missing(height)) height <- length(genes)/2
    if(!missing(filename)) pdf(filename,width = 20,height = height)

    if(missing(bottom_annotation)) {
        bottom_annotation <- NULL
    } else {
        annotation <- NULL
        bottom_annotation <-  HeatmapAnnotation(annotation_height = unit(0.3, "cm"),
                                                df = NULL,
                                                #col = list(Type = c("LGG"="green", "GBM"="orange")),
                                                annotation_legend_param = list(title_gp = gpar(fontsize = 16, fontface = "bold"),
                                                                               labels_gp = gpar(fontsize = 16), # size labels
                                                                               grid_height = unit(8, "mm")))
    }
    p <- oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                   row_order = NULL,
                   remove_empty_columns = FALSE,
                   column_order = NULL, # Do not sort the columns
                   alter_fun_list = alter_fun_list, col = color,
                   row_names_gp = gpar(fontsize = 16),  # set size for row names
                   pct_gp = gpar(fontsize = 16), # set size for percentage labels
                   axis_gp = gpar(fontsize = 16),# size of axis
                   #column_title = "OncoPrint for TCGA LGG, genes in Glioma signaling",
                   #column_title_gp = gpar(fontsize = 11),
                   row_barplot_width = unit(4, "cm"), #size barplot
                   bottom_annotation = bottom_annotation,
                   heatmap_legend_param = list(title = "Mutations", at = c("DEL", "INS", "SNP"),
                                               labels = c("DEL", "INS", "SNP"),
                                               title_gp = gpar(fontsize = 16, fontface = "bold"),
                                               labels_gp = gpar(fontsize = 16), # size labels
                                               grid_height = unit(8, "mm") # vertical distance labels
                   )
    )
    print(p)
    if(!missing(filename)) dev.off()
}
