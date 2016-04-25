![](inst/app/www/logo_red.png)

# Introduction

TCGAbiolinksGUI was created to help users without knowledge of programming to search, download and analyze 
TCGA data. This package offers an graphical user interface to the R/biocondcutor packages [TCGAbiolinks]( 	http://bioconductor.org/packages/TCGAbiolinks/)  and [ELMER](http://bioconductor.org/packages/ELMER/) packages.
Also, some other useful packages from bioconductor, such as [ComplexHeatmap](http://bioconductor.org/packages/ComplexHeatmap/)  package  has been used for data visualization.

## Installation TCGAbiolinksGUI

To install the package from biocondcutor repository, please, use the code below.

```R
source("https://bioconductor.org/biocLite.R")
biocLite("TCGAbiolinksGUI")
```

To install the package from a binary package, please, use the code below.

```R
# dependencies
devtools::install_github("BioinformaticsFMRP/TCGAbiolinks")
source("https://bioconductor.org/biocLite.R")
biocLite(c("pathview","clusterProfiler"))
install.packages(c("UpSetR","shiny","ReporteRs"))
devtools::install_github("thomasp85/shinyFiles")
devtools::install_github("ebailey78/shinyBS", ref="shinyBS3")
devtools::install_github("daattali/shinyjs")
install.packages("~/TCGAbiolinksGUI_0.99.0_R_x86_64-pc-linux-gnu.tar.gz", repos = NULL, type = "source")
```

## Quick start

The following commands should be used in order to start the graphical user interface.

```R
library(TCGAbiolinksGUI)
TCGAbiolinksGUI()
```

## Citation

Please cite both TCGAbiolinks package and TCGAbiolinksGUI: 

* Colaprico A, Silva TC, Olsen C, Garofano L, Cava C, Garolini D, Sabedot T, Malta TM, Pagnotta SM, Castiglioni I, Ceccarelli M, Bontempi G and Noushmehr H. "TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data." Nucleic acids research (2015): gkv1507.

Also, if you have used ELMER analysis please cite:

*  Yao, L., Shen, H., Laird, P. W., Farnham, P. J., & Berman, B. P. "Inferring regulatory element landscapes and transcription factor networks from cancer methylomes." Genome Biol 16 (2015): 105.
* Yao, Lijing, Benjamin P. Berman, and Peggy J. Farnham. "Demystifying the secret mission of enhancers: linking distal regulatory elements to target genes." Critical reviews in biochemistry and molecular biology 50.6 (2015): 550-573.
