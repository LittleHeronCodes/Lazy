
ent2sym <- function(genes) {
	if(!exists('geneInfo')) {
		geneInfo = read.table(paste0(DIR_RESOURCE, '/DB/gene_info.csv'), sep = ',', header=TRUE)
	}
	if(all(grepl('^[0-9]+$', genes))) {
		out = geneInfo$hgnc_symbol[match(genes, geneInfo$entrez)]		
	} else {
		out = geneInfo$entrez[match(genes, geneInfo$hgnc_symbol)]
	}
	return(out)
}

