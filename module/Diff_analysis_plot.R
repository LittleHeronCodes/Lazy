### Gene expression analysis plot pack ###

# differential expression results with up down tagged

resultsDEGTagging = function(resultDF, q.co, fc.co, q.cname='adj.P.Val', lfc.cname='logFC') {
	# resultDF: data frame with adjusted p value(adj.P.Val, logFC cname
	# if(!all(c('adj.P.Val', 'logFC', 'AveExpr') %in% colnames(results3))) {
	# 	stop('Column names do not include adj.P.Val, logFC, AveExpr.')
	# }

	# upIdx = which(resultDF$adj.P.Val < q.co & resultDF$logFC > log2(fc.co))
	# dnIdx = which(resultDF$adj.P.Val < q.co & resultDF$logFC < -log2(fc.co))
	# qvIdx = which(resultDF$adj.P.Val < q.co)
	upIdx = which(resultDF[,q.cname] < q.co & resultDF[,lfc.cname] >  log2(fc.co))
	dnIdx = which(resultDF[,q.cname] < q.co & resultDF[,lfc.cname] < -log2(fc.co))
	qvIdx = which(resultDF[,q.cname] < q.co)

	resultDF$dir = NA
	# resultDF$dir[qvIdx] = 'qc'
	resultDF$dir[upIdx] = 'up'
	resultDF$dir[dnIdx] = 'dn'

	return(resultDF)
}


## MA
plotMACustom = function(resultDF, title.prefix = '', q.co, fc.co) {
	if(!all(c('adj.P.Val', 'logFC', 'AveExpr', 'dir') %in% colnames(results3))) {
		stop('Column names do not include adj.P.Val, logFC, AveExpr, dir.')
	}

	upIdx = which(resultDF$dir == 'up')
	dnIdx = which(resultDF$dir == 'dn')
	qvIdx = which(resultDF$adj.P.Val < q.co)

	mtitle = sprintf('%s FC > %g, (up : %d, dn : %d)', title.prefix, fc.co, length(upIdx), length(dnIdx))
	plot(  logFC ~ AveExpr, data = results3, pch = 20, main = mtitle, col = 'grey80', cex = .5)
	points(logFC ~ AveExpr, data = resultDF[qvIdx,], pch = 20, col = 'grey50', cex = .5)
	points(logFC ~ AveExpr, data = resultDF[upIdx,], pch = 20, col = 'red', cex = .5)
	points(logFC ~ AveExpr, data = resultDF[dnIdx,], pch = 20, col = 'green', cex = .5)
	abline(h = log2(fc.co), col = 'firebrick')
	abline(h = -log2(fc.co), col = 'firebrick')
	# abline(v = avgCo, col = 'steelblue')

}

plotMACustom2 = function(resultDF, title.prefix = '', q.co, fc.co, q.cname='adj.P.Val', Ave.cname='AveExpr', lfc.cname='logFC') {
	if(!all(c('adj.P.Val', 'logFC', 'AveExpr', 'dir') %in% colnames(results3))) {
		stop('Column names do not include adj.P.Val, logFC, AveExpr, dir.')
	}

	upIdx = which(resultDF$dir == 'up')
	dnIdx = which(resultDF$dir == 'dn')
	qvIdx = which(resultDF[,q.cname] < q.co)

	A = log2(resultDF[,Ave.cname])
	M = resultDF[,lfc.cname]

	mtitle = sprintf('%s FC > %g, (up : %d, dn : %d)', title.prefix, fc.co, length(upIdx), length(dnIdx))
	plot(  A, M, pch = 20, main = mtitle, col = 'grey80', cex = .5)
	points(A[qvIdx], M[qvIdx], pch = 20, col = 'grey50', cex = .5)
	points(A[upIdx], M[upIdx], pch = 20, col = 'red', cex = .5)
	points(A[dnIdx], M[dnIdx], pch = 20, col = 'green', cex = .5)
	abline(h = log2(fc.co), col = 'firebrick')
	abline(h = -log2(fc.co), col = 'firebrick')
	# abline(v = avgCo, col = 'steelblue')

}


A = log2(resultDF$baseMean)
M = resultDF$log2FoldChange

upIdx = which(resultDF$log2FoldChange >= 1.5 & resultDF$padj < 0.01)
dnIdx = which(resultDF$log2FoldChange <=-1.5 & resultDF$padj < 0.01)
degIdx = which(abs(resultDF$log2FoldChange) >= 1.5 & resultDF$padj < 0.01)

mtitle = paste(mType, '- cont,', '(up :', length(upIdx), 'dn :', length(dnIdx), ')')

png(sprintf('fig/RNA-seq_MA/A549_%s_MA(pcg).png', mType))
plot(A,M,pch = 20, cex = .6, col = 'grey45', main = mtitle)
points(A[degIdx], M[degIdx], col = 'red', pch = 20, cex = .6)
abline(h=0, col = 'red')
# text(locator(2), c(paste('up :', length(upIdx)), paste('dn :', length(dnIdx))))
dev.off()

## Volcano


## PCA
## MDS



