## Gene Interaction
## Gene Interaction Analysis -- functions

library(reshape2)
library(parallel)
library(ggplot2)


# Gene-Target Matrix mean crispr scores
createGeneTargetMatrix <- function(mut2cell, crMat, ncore){
	cat('Creating gene-target matrix...\n')
	tcheck = proc.time()
	mut2gene=mclapply(mut2cell, function(cv) {
		mm = crMat[,which(colnames(crMat) %in% cv)]
		mv = apply(mm,1,mean,na.rm=TRUE)
		return(mv)
		}, mc.cores = ncore)
	maMat=do.call(cbind, mut2gene)
	print((proc.time() - tcheck)/60)
	print(dim(maMat))

	return(maMat)
}


# Main effect calculation modified function
mainEffectsMine <- function(x, queryNeg, targetNeg, eps=1e-04, maxiter=100) {

	cat('Calculating Main Effects...\n')
	tcheck = proc.time()
	# initialize
	z <- x
	nr <- nrow(z)
	nc <- ncol(z)
	r <- numeric(nr)
	c <- numeric(nc)

	neg <- 0
	oldsum <- 0

	# iteration
	for(i in 1:maxiter) {
		tcheck <- proc.time()

		rdelta <- apply(z, 1, median, na.rm = TRUE)	#**
		z <- z - matrix(rdelta, nrow=nr, ncol=nc, byrow=FALSE)
		r <- r + rdelta

		cdelta <- apply(z, 2, median, na.rm = TRUE)	#**
		z <- z - matrix(cdelta, nrow=nr, ncol=nc, byrow=TRUE )
		c <- c + cdelta

		if(missing(queryNeg)) {
			delc <- median(c, na.rm = TRUE)
		} else {
			delc <- mean(c[queryNeg], na.rm = TRUE)
		}
		if(!is.finite(delc)) delc = 0

		if(missing(targetNeg)) {
			delr <- median(r, na.rm = TRUE)
		} else {
			delr <- mean(r[targetNeg], na.rm = TRUE)
		}
		if(!is.finite(delr)) delr = 0

		c <- c - delc
		r <- r - delr

		neg <- neg + delc + delr

		newsum <- sum(abs(z), na.rm = TRUE)	#*
		converged <- newsum == 0 || abs(newsum - oldsum) < eps * newsum
        ptime <- (proc.time()-tcheck)[3]
		cat('iter:',i,'[neg:',neg,'delta r:',delr,'c:',delc,'err:',newsum,'runtime:',ptime,']\n')

		if(converged) break
		oldsum <- newsum
	}
	if (!converged) {
	    warning(gettextf("maineffects() did not converge in %d iterations", maxiter), domain = NA)
	}
	ans <- list(neg = neg, targetMainEffect = r, queryMainEffect = c, pi = z)
	print((proc.time() - tcheck)/60)
	
	return(ans)
}


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



