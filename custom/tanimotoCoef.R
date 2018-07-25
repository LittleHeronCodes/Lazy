tanimotoCoef = function(A, B) {
	# A, B is vector
	int = intersect(A, B)
	uni = union(A, B)
	return(length(int) / length(uni))
}

# enrichment factor function
enrichment <- function(x, y, total) {
	expected = length(x) * length(y) / length(total)
	observed = length(intersect(x, y))
	if(expected == 0) return(NA)
	return(observed / expected)
}
