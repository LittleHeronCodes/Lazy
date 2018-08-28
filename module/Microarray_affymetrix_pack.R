### AFFYMETRIX ARRAY DATA NORMALIZATION custom function pack ###

library(limma)
library(RMySQL)
library(simpleaffy)

# read and normalize to expression matrix
readCelAndNormalizeWrapper = function(gseAccno, dir.path='data/', celpath=NULL) {
	require(simpleaffy)

	# read in CEL files
	if(is.null(celpath)) celpath = sprintf('%s/%s/cel', dir.path, gseAccno)
	fns = list.celfiles(path = celpath, full.names = TRUE)
	print(fns)

	if(length(fns) == 0) stop('No CEL file detected')
	cat('Reading CEL files...\n')
	rawData = ReadAffy(filenames = fns)
	# define covdesc if covdesc file is saved under different name (no need)

	dataRMA = call.exprs(rawData, "rma")

	# extract expression matrix
	exprMat = exprs(dataRMA)
	colnames(exprMat) = gsub('.CEL', '', colnames(exprMat))

	return(exprMat)	
}


# Build model design matrix
buildDesignMatrix = function(model, name = 'group', levels=NULL) {
	groups = model[,name]
	if(is.null(levels)) levels = unique(groups)
	groups = factor(groups, levels = levels)
	design = model.matrix(~groups + 0)
	colnames(design) = levels(groups)

	return(design)
}


# Add psuedocounts (later)


# Linear fitting with limma package
lmFitEbayesWrapper = function(exprMat, design, contrMat, method = 'ls') {
	require(limma)

	fit = lmFit(exprMat, design, method = "ls")
	fit = eBayes(contrasts.fit(fit, contrMat))
	results = topTable(fit, number = nrow(exprMat), adjust = "BH", sort = "P")

	results$ID = rownames(results)
	results = results[,c(7,1:6)]
	return(results)
}


# Entrez gene mapping
entMapping = function(probeList, con, gpl) {
	require(RMySQL)
	pbs = paste0(probeList, collapse = "', '")
	stmt = paste("SELECT probe_id, entrez FROM probe2entrez",
		sprintf("WHERE gpl_id = '%s' AND probe_id IN ('%s');", gpl, pbs))
	genes = dbGetQuery(con, stmt)
	
	return(genes)
}


# DEG selection
probeSelect = function(results, q.co=.05, fc.co=2, avgCo=NULL) {
	# results : limma result data frame with probe level adj.P.Val, logFC, AveExpr
	if(is.null(avgCo)) avgCo = quantile(results2$AveExpr, 0.00)
	upIdx = which(results$adj.P.Val < q.co & results$logFC >  log2(fc.co) & results$AveExpr > avgCo)
	dnIdx = which(results$adj.P.Val < q.co & results$logFC < -log2(fc.co) & results$AveExpr > avgCo)
	print(paste('Probes up :', length(upIdx), 'dn :', length(dnIdx)))

	results$dir = NA
	results$dir[upIdx] = 'up'
	results$dir[dnIdx] = 'dn'

	return(results)
}

