![](inst/app/www/logo_gray2.png)

# Introduction

TCGAbiolinksGUI was created to help users without knowledge of programming to search, download and analyze 
TCGA data. This package offers an graphical user interface to the R/Bioconductor packages [TCGAbiolinks]( 	http://bioconductor.org/packages/TCGAbiolinks/)  and [ELMER](http://bioconductor.org/packages/ELMER/) packages.
Also, some other useful packages from Bioconductor, such as [ComplexHeatmap](http://bioconductor.org/packages/ComplexHeatmap/)  package  has been used for data visualization.

[Demo TCGAbiolinksGUI](https://tcgabiolinksgui.shinyapps.io/tcgabiolinks/)

## Installation TCGAbiolinksGUI

To install the package from [Bioconductor devel repository](http://bioconductor.org/packages/devel/bioc/html/TCGAbiolinksGUI.html), please, use the code below. The package will probably be in the Bioconductor release repository (stable) between April and May. 

```R
# for the moment it must be devel version of Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("TCGAbiolinksGUI", dependencies = TRUE)
```

To install the package development version from Github, please, use the code below.
```R
library(devtools)
source("https://bioconductor.org/biocLite.R", dependencies = T)

# dependencies
biocLite(c("GO.db", "DO.db","TCGAbiolinks"), dependencies = T)
devtools::install_github("BioinformaticsFMRP/TCGAbiolinksGUI", dependencies = T)
```

## Video tutorials

To facilitate the use of this package, we have created some tutorial videos demonstrating the tool.
Please check this [youtube list](https://www.youtube.com/playlist?list=PLoDzAKMJh15m40f7OqOLAW0nJwkVStJIJ).


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

