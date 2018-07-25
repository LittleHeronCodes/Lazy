## Calculate enrichment factor from data frame

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

	ef = ( cntPct[idx] * nrow(df) ) / ( length(rowIdx) * cntAll[idx] )
	ef = sort(ef, decreasing = TRUE)

	# add zero values
	idx2= setdiff(names(cntAll), idx)
	ef2 = rep(0, length(idx2))
	names(ef2) = idx2

	return(c(ef, ef2))
}
