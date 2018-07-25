# Validating Illumina raw data file format
## Adopted from CDRgator array_analysis

checkFormatIlmn <- function(fn) {
	ln <- system2("grep", c("'ILMN_'", fn, "-m 1 -n -o"), stdout = TRUE)
	if(length(ln) < 1) stop("Wrong or missing Probe Id. Check file.")
	sn <- as.numeric(substring(ln, 1,1))
	cm <- strsplit(system2("head", c(paste("-n",sn), fn), stdout = TRUE), '\t')

	id_column <- trimws(cm[[sn-1]][grep('ILMN_', cm[[sn]])])
	vch <- c(any(grepl('SAMPLE', cm[[sn-1]])), any(grepl('AVG_Signal', cm[[sn-1]])))
	vl_column <- c('SAMPLE', 'AVG_Signal')[vch]
	if(length(vl_column) < 1) stop("Column name for expression values should include 'ID_REF' or 'AVG_Signal")
	return(c(id_column, vl_column))
}
