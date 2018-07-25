# entrez gene mapping

entMapping <- function(probeList, con, platform) {
	require(RMySQL)
	pbs <- paste0(probeList, collapse = "', '")
	stmt <- sprintf("SELECT probe_id, entrez FROM probe2entrez WHERE gpl_id = '%s' AND probe_id IN ('%s');", platform, pbs)
	genes <- dbGetQuery(con, stmt)
	return(genes)
}

