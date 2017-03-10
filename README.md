![](inst/app/www/logo_gray2.png)

# Introduction

TCGAbiolinksGUI was created to help users without knowledge of programming to search, download and analyze 
TCGA data. This package offers an graphical user interface to the R/Bioconductor packages [TCGAbiolinks]( 	http://bioconductor.org/packages/TCGAbiolinks/)  and [ELMER](http://bioconductor.org/packages/ELMER/) packages.
Also, some other useful packages from Bioconductor, such as [ComplexHeatmap](http://bioconductor.org/packages/ComplexHeatmap/)  package  has been used for data visualization.

A running version of the GUI is found in (shinyapps.io)[https://tcgabiolinksgui.shinyapps.io/demo/]

## Installing TCGAbiolinksGUI

To install the package from [Bioconductor devel repository](http://bioconductor.org/packages/devel/bioc/html/TCGAbiolinksGUI.html), please, use the code below. The package will probably be in the Bioconductor release repository (stable) between April and May. 

```R
# for the moment it must be devel version of Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("TCGAbiolinksGUI", dependencies = TRUE)
```

To install the package development version from Github, please, use the code below.
```R
source("https://bioconductor.org/biocLite.R")
deps <- c("pathview","clusterProfiler","ELMER", "DO.db","GO.db", "ComplexHeatmap","EDASeq", "TCGAbiolinks")
for(pkg in deps)  if (!pkg %in% installed.packages()) biocLite(pkg, dependencies = TRUE)
deps <- c("devtools","shape","shiny","readr","googleVis","shinydashboard","shinyFiles","shinyjs","shinyBS")
for(pkg in deps)  if (!pkg %in% installed.packages())  install.packages(pkg,dependencies = TRUE)
devtools::install_github("BioinformaticsFMRP/TCGAbiolinksGUI")
```

## Docker image

TCGAbiolinksGUI is also available as Docker image (self-contained environments that contain everything needed to run the software), 
which can be easily run on Mac OS, Windows and Linux systems. 

This [PDF](https://drive.google.com/open?id=0B0-8N2fjttG-QXp5LVlPQnVQejg) show how to install and execute the image.

The image can be obtained from Docker Hub: https://hub.docker.com/r/tiagochst/tcgabiolinksgui/

For more information please check: https://docs.docker.com/ and https://www.bioconductor.org/help/docker/


## Video tutorials

To facilitate the use of this package, we have created some tutorial videos demonstrating the tool.
Please check this [youtube list](https://www.youtube.com/playlist?list=PLoDzAKMJh15m40f7OqOLAW0nJwkVStJIJ).

## PDF tutorials

For each section we created some PDFs with detailing the steps of each section: 
[Link to folder with PDFs](https://drive.google.com/drive/folders/0B0-8N2fjttG-Q25ldVVmUTVOTk0?usp=sharing)


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

