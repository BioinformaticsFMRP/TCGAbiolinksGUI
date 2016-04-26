#' @importFrom rvest html_table
#' @importFrom xml2 read_html
# @importFrom TCGAbiolinks TCGAquery
#' @keywords internal
load.maf <- function(env){

    tables <- xml2::read_html("https://wiki.nci.nih.gov/display/TCGA/TCGA+MAF+Files")
    tables <-  rvest::html_table(tables)
    # Table one is junk
    tables[[1]] <- NULL

    # get which tables are from the tumor
    all.df <- data.frame()
    for(tumor in unique(TCGAbiolinks::TCGAquery()$Disease)) {

        idx <- which(mapply(function(x) {
            any(grepl(tumor,(x[,1]), ignore.case = TRUE))
        },tables) == TRUE)
        df <- lapply(idx,function(x) tables[x])

        if(length(df) == 0) next
        # merge the data frame in the lists
        if(length(idx) > 1) {
            df <- Reduce(function(...) merge(..., all=TRUE), df)
        }  else if(length(idx) == 1) {
            df <- Reduce(function(...) merge(..., all=TRUE), df)
            df <- df[[1]]
            colnames(df) <- gsub(" ",".", colnames(df))
            colnames(df) <- gsub(":",".", colnames(df))
        }

        # Remove obsolete/protected
        df <- subset(df, df$Deploy.Status == "Available")
        df <- subset(df, df$Protection.Status == "Public")

        if(nrow(df) == 0) next

        df$Tumor <- tumor
        all.df <- rbind(all.df,df)
    }

    all.df[,"Deploy.Location"] <- gsub("/dccfiles_prod/tcgafiles/",
                                       "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/",
                                       all.df[,"Deploy.Location"] )
    assign("maf.files",all.df, envir = env)
}
