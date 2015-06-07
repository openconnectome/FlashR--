
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

dominate.greedy <- function (g)
{
    od <- degree(g, mode = "out") + 1
    S <- NULL
    A <- get.adjacency(g)
    diag(A) <- 0
    n <- nrow(A)
    covered <- rep(0, n)
    while (sum(covered) < n ) {
        i <- which.max(od)
        covered[A[i, ] > 0] <- 1
        covered[i] <- 1
        S <- c(S, i)
        A[, covered > 0] <- 0
        h <- graph.adjacency(A, mode = "directed")
        od <- degree(h, mode = "out") + 1 - covered
    }
    S
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

fastPlot <- function(g,layout,use.alpha,alpha)
{
   plot(layout,pch=20,axes=FALSE,xlab="",ylab="")
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
   points(layout,pch=20,col=2)
}

fastPlot3D <- function(g,layout)
{
   if(ncol(layout)==2) layout <- cbind(layout,rep(0,nrow(layout)))
   plot3d(layout,axes=FALSE,xlab="",ylab="",zlab="",box=FALSE)
   #edges <- get.edgelist(g,names=FALSE)
   #segments3d(layout[edges[,1],1],layout[edges[,1],2],layout[edges[,1],3],
   #          layout[edges[,2],1],layout[edges[,2],2],layout[edges[,2],3])
}

graph.spectral.embedding <- function(graph=g,no=2,
   cvec = degree(graph)/(vcount(graph) - 1))
{
	A <- get.adjacency(graph)
	B <- A+Diagonal(x=cvec)
	svds(B,k=no)
}

getLayout <- function(g,plotMethod,u, FRniter, FRcoolexp,
    circular, star.center,
   n, KKniter, KKinittemp, KKcoolexp)
{
  if(is.null(g)) return(NULL)
  layout <- paste('layout',gsub(" ",".",tolower(plotMethod)),
                  sep=".")
  if(layout == 'layout.auto'){
     layout <- layout.auto(g,dim=2)
  } else if(layout == 'layout.laplacian'){
     A <- graph.laplacian(as.undirected(g,mode='collapse'))
     z <- eigs(A,k=4,which="SM")
     d <- rev(z$values)
     layout <- z$vectors[,(length(d)-1):1]
  } else if(layout == 'layout.rdpg'){
     z <- graph.spectral.embedding(g,no=3)
     if(u=="U") {
        x <- z$u
     } else if(u=="V"){
        x <- z$v
     } else {
        x <- cbind(z$u[,1],z$v[,1])
     }
     layout <- x
  } else if(layout == 'layout.fruchterman.reingold'){
     layout <- layout.fruchterman.reingold(g,niter=FRniter,
                                          coolexp=FRcoolexp)
  } else if(layout == 'layout.fruchterman.reingold.grid'){
     layout <- layout.fruchterman.reingold.grid(g,niter=FRniter,
                                          coolexp=FRcoolexp)
  } else if(layout == 'layout.reingold.tilford'){
     layout <- layout.reingold.tilford(g,circular=circular)
  } else if(layout == 'layout.star'){
     layout <- layout.star(g,center=min(star.center,n))
  } else if(layout == 'layout.kamada.kawai'){
     layout <- layout.kamada.kawai(g,niter=KKniter,
                                          inittemp=KKinittemp,
                                          coolexp=KKcoolexp)
  } else if(layout =='layout.coordinates'){
     x <- get.vertex.attribute(g,'x')
     if(!is.null(x)){
        y <- get.vertex.attribute(g,'y')
        layout <- cbind(x,y)
     } else {
        layout <- layout.auto(g,dim=2)
     }
  } else {
     layout <- get(layout)(g)
  }
  layout

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

