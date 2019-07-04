# RNA sequencing count normalization functions

counts_to_TPM <- function(counts, featureLength) {

	stopifnot(length(featureLength) == nrow(counts))
	stopifnot(names(featureLength) == rownames(counts))

	x = counts / featureLength
	tpm = t( t(x) * 1e6 / colSums(x) )

	return(tpm)
}




