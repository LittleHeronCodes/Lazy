# Database connection function (server)

# connection function
perconn <- function(db = "", server = TRUE, usern = "", pw = "") {
	require(RMySQL)
	if(server) {
		conn <- dbConnect(dbDriver("MySQL"), host = "203.255.191.201", user = usern, password = pw, dbname = db, port = 3307)
	} else{
		conn <- dbConnect(dbDriver("MySQL"), user='', password='', dbname=db)
	}
	return(conn)
		}

# # for quick querying
# quickQuery <- function(con, stmt) {
# 	require(RMySQL)
# 	conn <- con
# 	results <- dbGetQuery(conn, stmt)
# 	dbDisconnect(con)
# 	return(results)
# }


