## make python style dictionary for indexing

makeDictionary = function(keys,values) {
	if(length(keys) != length(values)) stop('Keys and values must be of same length')
	idx = values
	names(idx) = keys
	return(idx)
}
