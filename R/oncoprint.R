#' Creating a oncoprint
#' @param mut Mutation file (see TCGAquery_maf from TCGAbiolinks)
#' @param genes Gene list
#' @param filename name of the pdf
#' @param color named vector for the plot
#' @param height pdf height
#' @importFrom ComplexHeatmap oncoPrint
#' @importFrom grid gpar grid.rect
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
                              height,
                              label.title = "Mutation",
                              font.size = 16){

    if(missing(mut))   stop("Missing mut argument")

    if(!missing(genes) & !is.null(genes)) mut <- subset(mut, mut$Hugo_Symbol %in% genes)
    mut$value <- 1
    mat <- dcast(mut, Tumor_Sample_Barcode + Hugo_Symbol ~ Variant_Type,value.var = "value",fill = 0)

    columns <- colnames(mat)[-c(1:2)]
    if(missing(color)){
        color <- c(rainbow(length(columns)))
        names(color) <- columns
    }
    mat$value <- ""
    for ( i in columns){
        mat[,i] <-  replace(mat[,i],mat[,i]>0,paste0(i,";"))
        mat[,i] <-  replace(mat[,i],mat[,i]==0,"")
        mat$value <- paste0(mat$value,mat[,i])
    }

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
    alt = intersect(names(alter_fun), c("background",as.character(columns)))
    alter_fun <- alter_fun[alt]
    mat <- t(mat)
    if(!missing(height)) height <- length(genes)/2
    if(!missing(filename)) pdf(filename,width = 20,height = height)

    p <- oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                   row_order = NULL,
                   remove_empty_columns = FALSE,
                   column_order = NULL, # Do not sort the columns
                   alter_fun = alter_fun, col = color,
                   row_names_gp = gpar(fontsize = font.size),  # set size for row names
                   pct_gp = gpar(fontsize = font.size), # set size for percentage labels
                   axis_gp = gpar(fontsize = font.size),# size of axis
                   #column_title = "OncoPrint for TCGA LGG, genes in Glioma signaling",
                   #column_title_gp = gpar(fontsize = 11),
                   row_barplot_width = unit(2, "cm"), #size barplot
                   #bottom_annotation = bottom_annotation,
                   heatmap_legend_param = list(title = label.title, at = names(color),
                                               labels = names(color),
                                               title_gp = gpar(fontsize = font.size, fontface = "bold"),
                                               labels_gp = gpar(fontsize = font.size), # size labels
                                               grid_height = unit(8, "mm") # vertical distance labels
                   )
    )
    print(p)
    if(!missing(filename)) dev.off()
}
