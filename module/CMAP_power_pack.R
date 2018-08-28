### CMAP WITH LINC CALCULATE TANIMOTO custom function pack ###

library(parallel)
library(pbmcapply)

## Required datasets
# LINC Signature : load(paste0(DIR_CMAP_Robj, '/trt_cp_total_siganture_top',lsig,'.Rdata'))
# sig_up, sig_dn, sinfo, pinfo, ginfo, sinfo, cinfo


# tanimoto coefficients
tanimotoCoef = function(A, B) {
	# A, B is vector
	int = intersect(A, B)
	uni = union(A, B)
	return(length(int) / length(uni))
}

# tanimoto similarity signature ranking
CMAPSaveTanimotoResultFiles <- function(geneList, lsig, odir, ncore) {

	if(any(!sapply(c('ginfo','sinfo','sig_up','sig_dn','tanimotoCoef'), exists)))
		stop('REQUIRED OBJECTS : ginfo, sinfo, sig_up, sig_dn, tanimotoCoef')

	print(sapply(geneList, function(ls) sapply(ls, length)))

	if(grepl('_bing', lsig)) {
		geneSpace = intersect(ginfo$pr_gene_id[which(ginfo$pr_is_bing == 1)], geneList[[1]]$toGene)
	} else geneSpace = intersect(ginfo$pr_gene_id, geneList[[1]]$toGene)

	cat('gene space:',length(geneSpace),'\n')

	if(!file.exists(odir)) {
		dir.create(odir, recursive = TRUE)
		cat('Directory',odir, 'created!\n')
	}

	for(mType in names(geneList)) {

		cat(mType, 'started!\n')
		upGene = intersect(geneList[[mType]]$upGene, geneSpace)
		dnGene = intersect(geneList[[mType]]$dnGene, geneSpace)

		# LINC signature selection
		tcheck = proc.time()
		tancoLs = pbmclapply(sinfo$sig_id, function(sigID) {

			# get linc signature
			lincSigLs.tmp = list(upLinc = intersect(sig_up[[sigID]], geneSpace),
								 dnLinc = intersect(sig_dn[[sigID]], geneSpace))

			# calculate tanimoto
			tancoUD = tanimotoCoef(upGene, lincSigLs.tmp$dnLinc)
			tancoDU = tanimotoCoef(dnGene, lincSigLs.tmp$upLinc)
			tancoUP = tanimotoCoef(upGene, lincSigLs.tmp$upLinc)
			tancoDN = tanimotoCoef(dnGene, lincSigLs.tmp$dnLinc)
			tancoM1 = (tancoUD + tancoDU) / 2
			tancoM2 = (tancoUP + tancoDN) / 2

			output.tmp = data.frame(
				sig_id = sigID,
				tanimoto.ud = tancoUD, tanimoto.du = tancoDU,
				tanimoto.uu = tancoUP, tanimoto.dd = tancoDN,
				tanimoto.mean1 = tancoM1, tanimoto.mean2 = tancoM2
				)

			return(output.tmp)
			}, mc.cores = ncore)

		names(tancoLs) = sinfo$sig_id
		print( (proc.time() - tcheck) / 60 )

		tcheck = proc.time()
		tanimotoResDF = do.call('rbind', tancoLs)
		print( (proc.time() - tcheck) / 60 )

		# save by cell line (file too big to save at once!)
		for(cell in unique(sinfo$cell_id)) {
			sigID.save = sinfo[which(sinfo$cell_id == cell),'sig_id']
			sinfo.save = sinfo[which(sinfo$sig_id %in% sigID.save),]
			tanimotoResDF.save = tanimotoResDF[which(tanimotoResDF$sig_id %in% sigID.save),]

			save(tanimotoResDF.save, sinfo.save, 
				file = sprintf('%s/%s_%ssig_%s_tanimotoResults.RData', odir, mType, lsig, cell))
			cat(mType, '-', cell,'saved!\n')
		}
	}
}


#==================================================================

# read tanimotoResults RData file
readTanimotoResultFiles <- function(mType, lsig, odir) {
	flist = list.files(odir, pattern = sprintf('^%s_%ssig_.*_tanimotoResults.RData', mType, lsig))
	tanimotoResDFLs = list()
	sinfoLs = list()
	for(fn in flist) {
		load(paste0(odir,'/', fn))
		tanimotoResDFLs[[fn]] = tanimotoResDF.save
		sinfoLs[[fn]] = sinfo.save
		cat(fn, 'loaded!\n')
	}
	tanimotoResDF = do.call('rbind', tanimotoResDFLs)
	sinfo = do.call('rbind', sinfoLs)

	rownames(tanimotoResDF) = NULL
	rownames(sinfo) = NULL
	return(list(resDF = tanimotoResDF, sinfo=sinfo)	)
}


## Calculate enrichment factor for drugs
enrichmentFactorForDataFrame = function(df, cutoff = 0.05, enrich_element, rank_by,psc=0) {
	##	df : data frame
	##	cutoff : cut percentile
	##	enrich_element : category column from df to enrich
	##	rank_by : column used to sort rank

	# order df by rank, get only within cut off
	co = sort(df[, rank_by], decreasing = TRUE)[nrow(df)*cutoff]
	rowIdx = which(df[,rank_by] >= co)
	# rowIdx = order(df[, rank_by], decreasing = TRUE)[1:(nrow(df)*cutoff)]

	# counts
	cntAll = table(df[, enrich_element])
	cntPct = table(df[rowIdx, enrich_element])
	idx = names(cntPct)

	ef = ( cntPct[idx] + psc ) / ( length(rowIdx) * cntAll[idx] / nrow(df) + psc)
	ef = sort(ef, decreasing = TRUE)
	rat = paste0(cntPct[idx], '/', cntAll[idx])
	names(rat) = idx
	rat = rat[names(ef)]

	# add zero values
	idx2= setdiff(names(cntAll), idx)
	ef2 = rep(0, length(idx2))
	names(ef2) = idx2
	rat2 = paste0(0, '/', cntAll[idx2])
	names(rat2) = idx2

	# rank = length(c(ef,ef2)) + 1 - rank(c(ef,ef2), ties.method='max')
	rank2 = rank(c(ef,ef2), ties.method='max') / length(c(ef,ef2)) * 100
	rank = 100/(100-min(rank2)) * (rank2 - min(rank2))

	return( list(ef = c(ef, ef2), ratio = c(rat,rat2), rank = rank) )
}


hypergeotestForDataFrame = function(df, cutoff = .05, enrich_element, rank_by) {
	##	df : data frame
	##	cutoff : cut percentile
	##	enrich_element : category column from df to enrich
	##	rank_by : column used to sort rank, column with rank percentile

	# order df by rank, get only within cut off
	rowIdx = order(df[, rank_by], decreasing = TRUE)[1:(nrow(df)*cutoff)]

	# counts
	U = nrow(df)


	cntAll = table(df[, enrich_element])
	cntPct = table(df[rowIdx, enrich_element])
	idx = names(cntPct)

	# idx = names(cntAll)

	g = idx[1]
	hg = phyper(cntPct[g], cntAll[g], nrow(df) - cntAll[g], length(rowIdx), lower.tail = TRUE)

	# add zero values
	idx2= setdiff(names(cntAll), idx)
	ef2 = rep(0, length(idx2))
	names(ef2) = idx2
	# ef = c(ef, ef2)

	return(c(ef, ef2))
}


## gene Count and change to symbol
geneCount = function(v){
	v = unlist(strsplit(v,split=";"))
	v = v[v!=""]
	if(length(v) == 0) return("")
	v = ginfo[match(v, ginfo[,1]),2]
	tb = table(v)
	tb = sort(tb, decreasing = TRUE)
	paste(paste(names(tb),"(",tb,")",sep=""),collapse=";")
}






