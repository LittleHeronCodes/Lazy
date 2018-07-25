### CMAP WITH LINC CALCULATE TANIMOTO custom function pack ###

library(parallel)
library(pbmcapply)
library(data.table)

## Required datasets
# LINC Signature : load(paste0(DIR_CMAP_Robj, '/trt_cp_total_siganture_top',lsig,'.Rdata'))
# sig_up, sig_dn, sinfo, pinfo, ginfo, sinfo, cinfo

# drug_summ_map = read.table(paste0(DIR_RESOURCE, '/DB/DATA__PHARMACOGX_DB_DRUGS_SUMMARY.csv'), sep = ',', header = TRUE)
# lincCompMap = drug_summ_map[,c(1,2,4,6,7,9,15)]
# lincCompMap = lincCompMap[which(lincCompMap$lincs_pert_id != ''),]


# tanimoto coefficients
tanimotoCoef = function(A, B) {
	# A, B is vector
	int = intersect(A, B)
	uni = union(A, B)
	return(length(int) / length(uni))
}


# calculating Tanimoto coef for each sig_id (INTERNAL FUNCTION ONLY)
CMAPGetTanimotoForSigID = function(sigID, upGene, dnGene, sig_up, sig_dn, intersectGenes=FALSE) {

	# get linc signature
	lincSigLs.tmp = list(upLinc = intersect(sig_up[[sigID]], geneSpace),
						 dnLinc = intersect(sig_dn[[sigID]], geneSpace))

	# calculate tanimoto
	tancoUD = tanimotoCoef(upGene, lincSigLs.tmp$dnLinc)
	tancoDU = tanimotoCoef(dnGene, lincSigLs.tmp$upLinc)
	tancoUp = tanimotoCoef(upGene, lincSigLs.tmp$upLinc)
	tancoDn = tanimotoCoef(dnGene, lincSigLs.tmp$dnLinc)
	tancoMean = (tancoUD + tancoDU) / 2

	if(intersectGenes) {
		# intersecting genes for target find
		intscUD = paste(intersect(upGene, lincSigLs.tmp$dnLinc), collapse = ';')
		intscDU = paste(intersect(dnGene, lincSigLs.tmp$upLinc), collapse = ';')		
	} else {
		intscUD = NA
		intscDU = NA
	}

	output.tmp = data.frame(
		sig_id = sigID,
		tanimoto.ud = tancoUD, tanimoto.du = tancoDU,
		tanimoto.uu = tancoUp, tanimoto.dd = tancoDn,
		tanimoto.mean = tancoMean,
		intersect.ud = intscUD,
		intersect.du = intscDU
		)

	return(output.tmp)

}

CMAPSaveTanCoResultsForGeneList = function(geneList, geneSpace, ncore, sinfo, odir, lsig, ...) {

	for(mType in names(geneList)) {

		cat(mType, 'started!\n')
		upGene = intersect(geneList[[mType]]$upGene, geneSpace)
		dnGene = intersect(geneList[[mType]]$dnGene, geneSpace)

		# LINC signature selection
		tcheck = proc.time()
		tancoLs = pbmclapply(sinfo$sig_id, function(sigID) {
			CMAPGetTanimotoForSigID(sigID, upGene=upGene,dnGene=dnGene, sig_up=sig_up, sig_dn=sig_dn)
			# CMAPGetTanimotoForSigID(sigID, upGene=upGene,dnGene=dnGene, ...)
			}, mc.cores = ncore)

		names(tancoLs) = sinfo$sig_id
		print( (proc.time() - tcheck) / 60 )

		tcheck = proc.time()
		tanimotoResDF = do.call('rbind', tancoLs)
		print( (proc.time() - tcheck) )

		tcheck = proc.time()
		tanimotoResDF = as.data.frame(rbindlist(tancoLs))
		print( (proc.time() - tcheck) )
		
		print(head(tanimotoResDF))

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
#==================================================================


## Calculate enrichment factor for drugs
enrichmentFactorForDataFrame = function(df, cutoff = 0.05, enrich_element, rank_by) {
	##	df : data frame
	##	cutoff : cut percentile
	##	enrich_element : category column from df to enrich
	##	rank_by : column used to sort rank

	# order df by rank, get only within cut off
	rowIdx = order(df[, rank_by], decreasing = TRUE)[1:(nrow(df)*cutoff)]

	# counts
	cntAll = table(df[, enrich_element])
	cntPct = table(df[rowIdx, enrich_element])
	idx = names(cntPct)

	ef = ( cntPct[idx] * nrow(df) + 1) / ( length(rowIdx) * cntAll[idx] + 1)
	ef = sort(ef, decreasing = TRUE)

	# add zero values
	idx2= setdiff(names(cntAll), idx)
	ef2 = rep(0, length(idx2))
	names(ef2) = idx2

	return(c(ef, ef2))
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






