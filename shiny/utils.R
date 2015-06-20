
html2txt <- function(str) {
		require(XML)
		if(nchar(str,type="bytes")==0) return(str)
		str <- paste("<html>",str,"</html>",sep="")
		paste(unlist(
			xpathApply(htmlParse(str, asText=TRUE),
						  "//body//text()", 
						  xmlValue)),
	      sep="\n",collapse="\n")
}

getOpenConnectomeList <- function()
{
   cat("Getting list of openconnectome graphs\n")

   openconnectome.graphs <- vector('list',length(openconnectome.animals))
   names(openconnectome.graphs) <- openconnectome.animals
   for(i in 1:length(openconnectome.graphs)){
      a <- scrape(paste(openconnectome.dir,openconnectome.animals[i],"/",sep=""),
                  parse=FALSE)
      x <- html2txt(a)
      b <- gregexpr("\n([[:alnum:]]|[\\._])+\\.graphml\\.?([[:alnum:]]+)?\n",x)
      openconnectome.graphs[[i]] <- gsub("\n","",unlist(regmatches(x,b)))
   }

   cat("Done\n")
   openconnectome.graphs

}


computeInvariants <- function(g,invariantList)
{
   cat("Computing Invariants:\n")
   inv <- data.frame(Invariants=invariants$Invariant,
                     Values=rep(0.0,nrow(invariants)),
                     Timings=rep(0.0,nrow(invariants)),
                     stringsAsFactors=FALSE)
   ind <- which(invariants$Invariant %in% invariantList)
   inv <- inv[ind,,drop=FALSE]
   invl <- invariants[ind,,drop=FALSE]
   for(i in 1:nrow(invl)){
      cat("Computing",invl$Function[i],"\n")
      t1 <- system.time(inv[i,2] <- 
          round(do.call(invl$Function[i],args=list(g=g))),3)
      inv[i,3] <- round(t1['elapsed'],3)
      cat("\t",inv[i,2],inv[i,3],"\n")
   }
   datatable(inv,rownames=FALSE)
}

fastPlot <- function(g,layout,use.alpha,alpha,size,color=2)
{
   plot(layout,pch=20,axes=FALSE,xlab="",ylab="",cex=size,col=color)
   edges <- get.edgelist(g,names=FALSE)
   if(use.alpha){
      col <- alpha('black',alpha)
   } else {
      col <- 1
   }
   if(is.directed(g)){
      arrows(layout[edges[,1],1],layout[edges[,1],2],
             layout[edges[,2],1],layout[edges[,2],2],length=0.1,col=col)
   } else {
      segments(layout[edges[,1],1],layout[edges[,1],2],
             layout[edges[,2],1],layout[edges[,2],2],col=col)
   }
   points(layout,pch=20,cex=size,col=color)
}

fastPlot3D <- function(g,layout,use.alpha,alpha,random)
{
   if(ncol(layout)==2) {
      if(random){
         layout <- cbind(layout,runif(nrow(layout)))
      } else {
         layout <- cbind(layout,rep(0,nrow(layout)))
      }
   }
   plot3d(layout,axes=FALSE,xlab="",ylab="",zlab="",box=FALSE,col=1)
   edges <- get.edgelist(g,names=FALSE)
   x <- matrix(0,nrow=2*nrow(edges),ncol=3)
   x[(1:nrow(edges))*2-1,] <- layout[edges[,1],]
   x[(1:nrow(edges))*2,] <- layout[edges[,2],]
   cat("plotting edges\n")
   if(use.alpha){
      segments3d(x[,1],x[,2],x[,3],col='black',alpha=alpha)
   } else {
      segments3d(x[,1],x[,2],x[,3],col=1)
   }
   cat("Done\n")
}

graph.spectral.embedding <- function(graph=g,no=2,
   cvec = degree(graph)/(vcount(graph) - 1))
{
	A <- get.adjacency(graph)
	B <- A+Diagonal(x=cvec)
	svds(B,k=no)
}

computeLaplacian <- function(graph,d,normalize)
{
     A <- graph.laplacian(as.undirected(graph,mode='collapse'),
              normalized=normalize)
     z <- eigs(A,k=d+1,which="SM")
     z$vectors[,1:d]
}

## code for invariant calculations
## these must take a graph and return a single number

max.component.size <- function(g){
   max(clusters(g)$csize)
}

get.girth <- function(g) girth(g)$girth
get.cent <- function(g) centralization.degree(g)$centralization
get.rec <- function(g){
   if(is.directed(g)) return(reciprocity(g))
   else return(NA)
}
get.fastgreedy <- function(g) length(fastgreedy.community(
            as.undirected(simplify(g))))

get.pl <- function(g) power.law.fit(
             degree.distribution(g))$alpha
get.core <- function(g) max(graph.coreness(g))
get.avg.degree <- function(g) mean(degree(g))
get.max.degree <- function(g) max(degree(g))
get.min.degree <- function(g) min(degree(g))
get.var.degree <- function(g) var(degree(g))
dominate.num.greedy <- function(g) length(dominate.greedy(g))

