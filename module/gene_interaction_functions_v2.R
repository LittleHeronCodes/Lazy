## Gene Interaction
## Gene Interaction Analysis -- functions

library(reshape2)
library(parallel)
library(ggplot2)

##=============================##
##  Analysis function wrapper  ##
##=============================##

# Gene-Target Matrix mean crispr scores (removed)

# Main effect calculation modified function (removed)


# t-test p values (NA exists)
getTestPValues2 <- function(gfMat, crMat, ncore){
	cat('Calculating p values...\n')
	tcheck = proc.time()
	pvMat = mclapply(rownames(gfMat), function(g) {
		midx = gfMat[which(rownames(gfMat) == g),] == 1
		ex1 = names(which(apply(crMat[,which(midx)],  1, function(v) sum(!is.na(v))) <= 2))
		ex2 = names(which(apply(crMat[,which(!midx)], 1, function(v) sum(!is.na(v))) <= 2))
		crMat.f = crMat[which(!rownames(crMat) %in% union(ex1, ex2)),]
		tres = apply(crMat.f, 1, function(v) t.test(v[which(midx)], v[which(!midx)], alternative='two.sided')$p.value)
		}, mc.cores = ncore)
	names(pvMat) = rownames(gfMat)
	pvMat = do.call(cbind, pvMat)
	print((proc.time() - tcheck)/60)
	# print(dim(pvMat))

	return(pvMat)
}

# t-test p values (NA doesn't exist)
getTestPValues <- function(gfMat, crMat, ncore){
	cat('Calculating p values...\n')
	tcheck = proc.time()
	pvMat = mclapply(rownames(gfMat), function(g) {
		midx = gfMat[which(rownames(gfMat) == g),] == 1
		tres = apply(crMat, 1, function(v) t.test(v[which(midx)], v[which(!midx)], alternative='two.sided')$p.value)
		}, mc.cores = ncore)
	names(pvMat) = rownames(gfMat)
	pvMat = do.call(cbind, pvMat)
	print((proc.time() - tcheck)/60)
	# print(dim(pvMat))

	return(pvMat)
}


# median difference (NA doesn't exist)
getMedDiff <- function(gfMat, crMat, ncore){
	cat('Calculating p values...\n')
	tcheck = proc.time()
	mdMat = mclapply(rownames(gfMat), function(g) {
		i = which(rownames(gfMat) == g)
		med1 = apply(crMat[,which(gfMat[g,] == 1)], 1, median)
		med2 = apply(crMat[,which(gfMat[g,] == 0)], 1, median)
		medD = med1-med2
		}, mc.cores = ncore)
	names(mdMat) = rownames(gfMat)
	mdMat = do.call(cbind, mdMat)
	print((proc.time() - tcheck)/60)
	# print(dim(mdMat))

	return(mdMat)
}


# mean difference (NA doesn't exist)
getMeanDiff <- function(gfMat, crMat, ncore){
	cat('Calculating p values...\n')
	tcheck = proc.time()
	mnMat = mclapply(rownames(gfMat), function(g) {
		i = which(rownames(gfMat) == g)
		med1 = apply(crMat[,which(gfMat[g,] == 1)], 1, mean)
		med2 = apply(crMat[,which(gfMat[g,] == 0)], 1, mean)
		medD = med1-med2
		}, mc.cores = ncore)
	names(mnMat) = rownames(gfMat)
	mnMat = do.call(cbind, mnMat)
	print((proc.time() - tcheck)/60)
	# print(dim(mnMat))

	return(mnMat)
}


# t-test p values (NA doesn't exist)
getKSTestPValues <- function(gfMat, crMat, ncore){
	cat('Calculating p values...\n')
	tcheck = proc.time()
	pvMat = mclapply(rownames(gfMat), function(g) {
		midx = gfMat[which(rownames(gfMat) == g),] == 1
		tres = apply(crMat, 1, function(v) ks.test(v[which(midx)], v[which(!midx)], alternative='two.sided')$p.value)
		}, mc.cores = ncore)
	names(pvMat) = rownames(gfMat)
	pvMat = do.call(cbind, pvMat)
	print((proc.time() - tcheck)/60)

	return(pvMat)
}







# t-test p vaules and median for drug-gene testing (NA exists)
tTestDrugGene <- function(gfMat, drMat, ncore){
	cat('Calculating p values...\n')
	tcheck = proc.time()
	outLs = mclapply(rownames(gfMat), function(g) {
		midx = gfMat[which(rownames(gfMat) == g),] == 1
		ex1 = names(which(apply(drMat[,which(midx)],  1, function(v) sum(!is.na(v))) <= 2))
		ex2 = names(which(apply(drMat[,which(!midx)], 1, function(v) sum(!is.na(v))) <= 2))
		drMat.f = drMat[which(!rownames(drMat) %in% union(ex1, ex2)),]

		tObj = apply(drMat.f, 1, function(v) t.test(v[which(midx)], v[which(!midx)], alternative='two.sided'))
		
		tpval = sapply(tObj, function(ls) ls$p.value)
		tstat = sapply(tObj, function(ls) ls$statistic)
		names(tstat) = gsub('.t$','', names(tstat))
		medDiff = apply(drMat.f, 1, function(v) median(v[which(midx)], na.rm=T) - median(v[which(!midx)], na.rm=T) )

		exfill = rep(NA, length(union(ex1,ex2)))
		names(exfill) = union(ex1,ex2)
		tstat = c(tstat, exfill)
		tstat = tstat[rownames(drMat)]
		tpval = c(tpval, exfill)
		tpval = tpval[rownames(drMat)]
		medDiff = c(medDiff, exfill)
		medDiff = medDiff[rownames(drMat)]

		return(list(tpval = tpval, tstat = tstat, medDiff = medDiff) )
		}, mc.cores = ncore)
	names(outLs) = rownames(gfMat)
	print((proc.time() - tcheck)/60)

	cat('Creating matrices...\n')
	drPvMat = do.call(cbind, lapply(outLs, function(ls) ls$tpval) )
	drMDMat = do.call(cbind, lapply(outLs, function(ls) ls$medDiff) )
	drTSMat = do.call(cbind, lapply(outLs, function(ls) ls$tstat) )

	# cat('Creating data frame...\n')
	# drAltDF = dgiResLs2DataFrame(outLs)

	return(list(drPvMat = drPvMat, drTSMat = drTSMat, drMDMat = drMDMat, drAltLs=outLs))
}

# change dgiRes list output from tTestDrugGene to easier access dataframe
dgiResLs2DataFrame <- function(outLs) {
	dfls = lapply(names(outLs), function(g) {
			tdf = data.frame(outLs[[g]])
			tdf$chem = rownames(tdf)
			tdf$altGene = g
			return(tdf)	})
	drAltDF = do.call(rbind, dfls)
	return(drAltDF)	
}


##======================##
##  Plotting Functions  ##
##======================##

# Volcano plots with known interaction

makePlotForKnownInt <- function(goi, known, piMat, pvMat, qvMat, draw='p') {

	# print(dim(piMat))
	if(missing(qvMat)) qvMat = apply(pvMat, 2, function(v) p.adjust(v, method = 'fdr') )

	piv = piMat[,which(colnames(piMat) == goi)]
	pvv = pvMat[,which(colnames(pvMat) == goi)]
	qvv = qvMat[,which(colnames(qvMat) == goi)]

	if(!identical(names(piv), names(pvv))) stop('Not equal')

	plotDF = data.frame( target_sym = names(piv), piScore = piv, t.testp = pvv, 
		logP = -log10(pvv), t.testq = qvv, logQ = -log10(qvv))

	plotDF$col = 'grey'
	plotDF$col[which(plotDF$target_sym %in% knownint)] = 'blue'
	plotDF$label = plotDF$target_sym
	plotDF$label[which(plotDF$col =='grey')] = ''

	if(draw=='p') {
		p = ggplot(plotDF, aes(x=logP, y=piScore))
	} else if(draw=='q') {
		p = ggplot(plotDF, aes(x=logQ, y=piScore))
	}

	p + 
	  geom_point(data=subset(plotDF, col == 'grey'),colour='grey') +
	  geom_point(data=subset(plotDF, col != 'grey'),colour='blue')	+
	  ggtitle(paste(goi, draw, 'value'))

}

# known interaction labelled
# knownint = giscores$Target_gene[which(giscores$Query_gene == goi & giscores$FDR_IHW < .1)]
# makePlotForKnownInt(goi, knownint, piMat, pvMat, qvMat,draw='q') +
#   geom_text(aes(label=label),vjust = 0, hjust=0, colour='blue')



