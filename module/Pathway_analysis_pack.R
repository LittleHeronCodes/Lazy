### PATHWAY ANALYSIS custom function pack ###

.libPaths(c(.libPaths(), RLIB))

library(clusterProfiler)
library(org.Hs.eg.db)
library(GSEABase)
library(RColorBrewer)

# available key types
print(keytypes(org.Hs.eg.db))

# read data from msigdb
c2.kegg = read.gmt(paste0(DIR_MSIG, '/c2.cp.kegg.v6.1.entrez.gmt'))
c2.bioc = read.gmt(paste0(DIR_MSIG, '/c2.cp.biocarta.v6.1.entrez.gmt'))
c5.gobp = read.gmt(paste0(DIR_MSIG, '/c5.bp.v6.1.entrez.gmt'))
c5.gocc = read.gmt(paste0(DIR_MSIG, '/c5.cc.v6.1.entrez.gmt'))
h.all   = read.gmt(paste0(DIR_MSIG, '/h.all.v6.1.entrez.gmt'))
pathways_db = c('c2.kegg', 'c2.bioc', 'c5.gobp', 'c5.gocc', 'h.all')
cat('\nLoaded gmts:', pathways_db, '\n')

#=====================================================
# pathway enrichment analysis
enricherForGeneListWrapper = function(glist, term, pcut=.1, qcut=.2) {
	# glist: gene list for one set with upGene, dnGene, toGene
	require(clusterProfiler)

	gspace = glist$toGene
	enrGeneLs = list(
		up   = glist$upGene,
		down = glist$dnGene,
		both = unique(unlist(glist[c('upGene', 'dnGene')]))
	)

	enr_res = lapply(enrGeneLs, function(gset) {
		output = enricher(gset, TERM2GENE = base::get(term), universe = gspace, 
			pvalueCutoff = pcut, pAdjustMethod = 'fdr', qvalueCutoff = qcut)
		return(output)
		})
	
	print(sapply(enr_res, nrow))
	return(enr_res)
}



# setting colors for heatmap by direction
pathwayAnalysisColorSet = function(dir) {
	require(RColorBrewer)
	if(!dir %in% c('up','down','both')) stop('dir needs to be one of up, down, both.')
	if(dir == 'up'	) colors = colorRampPalette(c('lemonchiffon','firebrick'))(101)
	if(dir == 'down') colors = colorRampPalette(c('lemonchiffon','seagreen4'))(101)
	if(dir == 'both') colors = colorRampPalette(c('azure','navy'))(101)
	return(colors)
}

pathwayAnalysisColorSet2 = function(dir) {
	require(RColorBrewer)
	if(!dir %in% c('up','down','both')) stop('dir needs to be one of up, down, both.')
	if(dir == 'up'	) colors = colorRampPalette(c('ivory','#890000'))(299)
	if(dir == 'down') colors = colorRampPalette(c('ivory','#003f01'))(299)
	if(dir == 'both') colors = colorRampPalette(c('azure','navy'))(299)
	return(colors)
}

pathwayAnalysisColorSet3 = function(dir) {
	require(circlize)
	# for ComplexHeatmap
	if(!dir %in% c('up','down','both')) stop('dir needs to be one of up, down, both.')
	if(dir == 'up'	) colors = colorRamp2(c(0,10), c('ivory','#890000'))
	if(dir == 'down') colors = colorRamp2(c(0,10), c('ivory','#003f01'))
	if(dir == 'both') colors = colorRamp2(c(0,10), c('azure','navy'))
	return(colors)
}

