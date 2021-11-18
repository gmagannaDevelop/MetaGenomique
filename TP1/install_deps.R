# Install bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", quietly = TRUE)
BiocManager::install("qvalue")

# Install pacakges
list.of.packages <-
  c('devtools', 'RJSONIO','ecodist','gplots', 'scatterplot3d','usethis', 'httr',
                  'rcmdcheck', 'roxygen2', 'rversions')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Install from github
library(devtools)
install_github(repo='MG-RAST/matR',quiet=T)
