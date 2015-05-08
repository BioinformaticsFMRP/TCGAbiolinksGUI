#' biOmics
#'
#' biOmics allows you to Download data of samples from TCGA, ENCODE, ROADMAP,
#' using a ontology-based algorithm
#'
#' The functions you're likely to need from \pkg{biOmics} are
#' \code{\link{biOmicsApp}}, \code{\link{geoDownloaderLinks}} and
#' \code{\link{filterRmapData}}. Otherwise refer to the vignettes to see
#' how to format the documentation.
#'
#' @docType package
#' @name biOmics
NULL

#' Terms of encode mapped automatically using ontology algorithm,
#' plus manual mapping for terms not found
#' @docType data
#' @keywords internal
#' @name biosample.encode
#' @format A data frame with 374 rows and 3 variables
NULL

#' Terms of tcga mapped automatically using ontology algorithm,
#' plus manual mapping for terms not found
#' @docType data
#' @keywords internal
#' @name biosample.tcga
#' @format A data frame with 35 rows and 3 variables
NULL

#' Terms of roadmap mapped automatically using ontology algorithm,
#' plus manual mapping for terms not found
#' @docType data
#' @keywords internal
#' @name biosample.roadmap
#' @format A data frame with 437 rows and 3 variables
NULL

#' TCGA disease table
#' @docType data
#' @keywords internal
#' @name disease.table
#' @format A data frame with 37 rows and 4 variables
NULL

#' TCGA platforms table
#' @docType data
#' @keywords internal
#' @name platform.table
#' @format A data frame with 79 rows and 4 variables
NULL

#' Platforms of roadmap, encode and TCGA mapped into groups.
#' @docType data
#' @keywords internal
#' @name platforms
#' @format A data frame with 123 rows and 2 variables
NULL

#' BTO values mapped for some selected systems.
#' @docType data
#' @keywords internal
#' @name systems
#' @format A data frame with 17 rows and 2 variables
NULL

#' Roadmap metadata database.
#' The data set contains the following fields:
#' \itemize{
#' \item accession Acession name
#' \item  biosample biosample name
#' \item assay Assay name
#' \item lab lab name
#' \item target target name
#' \item description description name
#' \item organism organim name
#' }
#' @docType data
#' @keywords internal
#' @name encode.db
#' @format A data frame with 4598 rows and 7 variables
NULL

#' Roadmap metadata database.
#' The data set contains the following fields:
#' \itemize{
#'   \item X..GEO.Accession Geo Acession number
#'   \item Sample.Name	Sample name
#'   \item Experiment	 Experiment name
#'   \item NA.Accession	NA acession number
#'   \item Center	Center name
#'   \item SRA.FTP	SRA file FTP path
#'   \item GEO.FTP	GEO files FTO path
#'   \item Embargo.end.date Embargo date
#'}
#' @docType data
#' @keywords internal
#' @name roadmap.db
#' @format A data frame with 3135 rows and 8 variables
NULL

#' TCGA metadata database.
#' The data set contains the following fields:
#' \itemize{
#' \item addedDate Date when sample was added to database
#' \item baseName	name of the sample folder
#' \item deployLocation	Path of the sample folder
#' \item deployStatus	If available or not
#' \item id TCGA id
#' \item isLatest Is the lasted version of the sample?
#' \item name Sample name
#' \item revision Sample revision
#' \item serialIndex Sample serial Index
#' \item Platform Sample platform
#' \item Disease Sample disease
#'}
#' @docType data
#' @keywords internal
#' @name tcga.db
#' @format A data frame with 4733 rows and 11 variables
NULL
