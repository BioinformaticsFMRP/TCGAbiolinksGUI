![](inst/app/www/logo_gray2.png)

# Introduction

TCGAbiolinksGUI was created to help users without knowledge of programming to search, download and analyze 
TCGA data. This package offers an graphical user interface to the R/biocondcutor packages [TCGAbiolinks]( 	http://bioconductor.org/packages/TCGAbiolinks/)  and [ELMER](http://bioconductor.org/packages/ELMER/) packages.
Also, some other useful packages from bioconductor, such as [ComplexHeatmap](http://bioconductor.org/packages/ComplexHeatmap/)  package  has been used for data visualization.

[Demo TCGAbiolinksGUI](https://tcgabiolinksgui.shinyapps.io/tcgabiolinks/)

## Installation TCGAbiolinksGUI

To install the package from biocondcutor repository, please, use the code below.

```R
source("https://bioconductor.org/biocLite.R")
biocLite("TCGAbiolinksGUI")
```

To install the package development version from Github, please, use the code below.
```R
library(devtools)
source("https://bioconductor.org/biocLite.R")

# dependencies
install.packages(c("shiny","readr","googleVis","shinydashboard"))
biocLite(c("pathview","clusterProfiler","ELMER", "GO.db", "DO.db"))
devtools::install_github("BioinformaticsFMRP/TCGAbiolinks")
devtools::install_github("thomasp85/shinyFiles")
devtools::install_github("ebailey78/shinyBS", ref="shinyBS3")
devtools::install_github("daattali/shinyjs")
devtools::install_github("BioinformaticsFMRP/TCGAbiolinksGUI")
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

* TCGAbiolinksGUI: A Graphical User Interface to analyze TCGA data. Manuscript in preparation.

Also, if you have used ELMER analysis please cite:

* Yao, L., Shen, H., Laird, P. W., Farnham, P. J., & Berman, B. P. "Inferring regulatory element landscapes and transcription factor networks from cancer methylomes." Genome Biol 16 (2015): 105.
* Yao, Lijing, Benjamin P. Berman, and Peggy J. Farnham. "Demystifying the secret mission of enhancers: linking distal regulatory elements to target genes." Critical reviews in biochemistry and molecular biology 50.6 (2015): 550-573.


If you have used  OncoPrint plot and Heatmap Plot please cite:

* Gu, Zuguang, Roland Eils, and Matthias Schlesner. "Complex heatmaps reveal patterns and correlations in multidimensional genomic data." Bioinformatics (2016): btw313

If you have used  Pathway plot please cite:

* Luo, Weijun, Brouwer and Cory (2013). “Pathview: an R/Bioconductor package for pathway-based data integration and visualization.” Bioinformatics, 29(14), pp. 1830-1831.

