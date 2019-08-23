# Database connection function (server)

# connection function
# perconn2 <- function(db = "", usern = "", pw = "", host = "", port=22) {
# 	require(RMySQL)
# 	conn = dbConnect(dbDriver("MySQL"), host = host, user = usern, password = pw, dbname = db, port = port)
# 	return(conn)
# }

# for quick querying
quickQuery <- function(con, stmt) {
	require(RMySQL)
	conn = con
	results = dbGetQuery(conn, stmt)
	dbDisconnect(con)
	return(results)
}

dbMatchList <- function(con, table, select_col, match_col, match_list) {
	select_col = paste0(select_col, collapse = ", ")
	match_list = paste0(match_list, collapse = "', '")

	stmt = sprintf("SELECT %s FROM %s WHERE %s IN ('%s');",	select_col, table, match_col, match_list)
	# print(stmt)
	results = dbGetQuery(con, stmt)
	return(results)
}


