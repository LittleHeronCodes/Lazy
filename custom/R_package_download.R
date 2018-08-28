## important packages in case of new installations

# install via cran
install.packages('BH')
install.packages('RMySQL')
install.packages('RColorBrewer')
install.packages('ggplot2')
install.packages('data.table')
install.packages('devtools')
install.packages('vennDiagram')

# tidyverse series
# install.packages('tidyverse')
install.packages('dplyr')
install.packages('reshape')
install.packages('magrittr')
install.packages('readxl')


# install via bioconductor
source('https://bioconductor.org/biocLite.R')

biocLite()
biocLite('BioInstaller')
biocLite('Biobase')
biocLite(c('Biostrings', 'BiocGenerics', 'BiocParallel'))
biocLite('DBI')
biocLite('AnnotationDbi')
biocLite('affy')
biocLite('gcrma')
biocLite('genefilter')
biocLite('simpleaffy')
biocLite('limma')
biocLite(c('hgu133plus2cdf', 'hgu133acdf', 'hgu133a2cdf','hugene10stv1cdf'))
biocLite('DESeq2')

biocLite('ComplexHeatmap')
biocLite(c("RBGL","graph")) # for Vennerable package


# install via devtools
library(devtools)
install_github("js229/Vennerable")
