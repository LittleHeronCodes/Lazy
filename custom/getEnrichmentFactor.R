# enrichment factor function
getEnrichmentFactor <- function(setA,setB,setT) {
	ef = NA
	A = length(setA)
	B = length(setB)
	T = length(setT)
	I = length(intersect(setA, setB))
	if((A / T * B) != 0) ef = I / (A / T * B)
	return(ef)
}
