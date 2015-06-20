

getOpenConnectome <- function(openconnectome.dir)
{
   a <- try(scrape(openconnectome.dir),silent=TRUE) 
   if(inherits(a,'try-error')){
      cat("Could not access openconnecto.me. Using cached data.")
      load("openConnectomeGraphs.RData")
      return(openconnectome.graphs)
   }
   ## get the first level of directories -- the species
   dirs <- xpathSApply(a[[1]],"//table//td/a",xmlValue)
   dirs <- dirs[grep('Parent Directory',dirs,invert=TRUE)]
   tree <- vector('list',length(dirs))
   h <- which(dirs=='human/')
   dirs <- c(dirs[-h],'human/')
   names(tree) <- gsub("/","",dirs)
   for(i in 1:length(dirs)){
      dir <- dirs[i]
      b <- scrape(paste(openconnectome.dir,dir,sep="/"))
      graphs <- xpathSApply(b[[1]],"//table//td/a",xmlValue)
      graphs <- graphs[grep('Parent Directory',graphs,invert=TRUE)]
      tree[[i]] <- graphs
   }
   tree
}
