# RNA sequencing count normalization functions

counts_to_TPM <- function(counts, featureLength) {

	stopifnot(length(featureLength) == nrow(counts))
	stopifnot(names(featureLength) == rownames(counts))

	tpm <- apply(counts, 2, function(x) {
		rpk = x / featureLength
		scalingFactor = sum(rpk) / 1e6
		return(rpk / scalingFactor)
	})

	return(tpm)
}

