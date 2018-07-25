
source('https://bioconductor.org/biocLite.R')

install.packages('BH')
install.packages('dplyr')
install.packages('reshape')
install.packages('RMySQL')
install.packages('RColorBrewer')
install.packages('ggplot2')

biocLite()
biocLite('AnnotationDbi')
biocLite('DBI')
biocLite('Biobase')
biocLite(c('Biostrings', 'BiocGenerics', 'BiocParallel'))
biocLite('affy')
biocLite('gcrma')
biocLite('genefilter')
biocLite('simpleaffy')
biocLite(c('hgu133plus2cdf', 'hgu133acdf', 'hgu133a2cdf','hugene10stv1cdf'))
biocLite('limma')


biocLite('BioInstaller')
