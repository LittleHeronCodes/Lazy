## Base helper functions

setcolfirst = function(DT, ...){
	nm = as.character(substitute(c(...)))[-1L]
	setcolorder(DT, c(nm, setdiff(names(DT), nm)))
}
