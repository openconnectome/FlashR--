rdf <- function(dist, ids, scans=2) {
  N <- dim(dist)[1]
  if (is.null((N))) {
    stop('Invalid datatype for N')
  }
  rdf <- array(NaN, N*(scans-1))
  count <- 1
  for (i in 1:N) {
    ind <- which(ids==ids[i])
    for (j in ind) {
      if (j != i) {
        di <- dist[i,]
        di[ind] <- Inf
        d <- dist[i,j]
        rdf[count] <- 1 - (sum(di[!is.nan(di)] < d) + 0.5*sum(di[!is.nan(di)] == d)) / (N-length(ind))
        count <-  count + 1
      }
    }
  }
  return(rdf)
}

mnr <- function(rdfs){
  mnr <- mean(rdfs[!is.nan(rdfs)])
}
