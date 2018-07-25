# fishers method sum p values

fishersMethod = function(x) pchisq(-2 * sum(log(x)),df = 2*length(x),lower=FALSE)

