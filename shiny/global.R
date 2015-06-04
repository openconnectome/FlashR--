library(igraph)
library(shinyBS)
library(scrapeR)
library(rARPACK)
library(Matrix)
library(DT)
library(scales)
library(FlashGraphR)


options(shiny.maxRequestSize = 300*1024^2)

plot.methods <- c('Auto','Random','circle',
                  'sphere','Fruchterman Reingold',
                  'Fruchterman Reingold Grid',
                  'Kamada Kawai',
                  'Spring',
                  'Reingold Tilford',
                  'LGL',
                  'star',
                  'Graphopt',
                  'RDPG',
                  'Laplacian',
                  'Coordinates')

plot.invariants <- c("Degree Distribution",
                "Alpha Centrality",
                "Betweenness",
                "Closeness",
                "Eigenvalue Centrality",
                "Boncich Power Centrality",
                "Authority Score",
                "Hub Score",
                "Articulation Points",
                "Burt's Constraint",
                "Diversity",
                "Average KNN Degree",
                "Page Rank")


openconnectome.dir <- "http://openconnecto.me/data/public/graphs/"
openconnectome.animals <- c("cat","fly","macaque","mouse","rat",
                            "worm","human")

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

invariants <- data.frame(Invariant=c(
                           "Order",
                           "Size",
                           "#Components",
                           "Max Component Size",
                           "Diameter",
                           "Girth",
                           "Density",
                           "Average Path Length",
                           "Centralization Degree",
                           "Reciprocity",
                           "Transitivity",
                           "Degree Assortativity",
                           "Maximal Coreness",
                           "Fast Greedy Communities",
                           "Power Law Fit (Exponent)"),
                        Function=c(
                           'vcount',
                           'ecount',
                           'no.clusters',
                           'max.component.size',
                           'diameter',
                           'get.girth',
                           'graph.density',
                           'average.path.length',
                           'get.cent',
                           'get.rec',
                           'transitivity',
                           'get.core',
                           'assortativity.degree',
                           'get.fastgreedy',
                           'get.pl'),
                       stringsAsFactors=FALSE)
                           

                         
