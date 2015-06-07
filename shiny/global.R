library(igraph)
library(shinyBS)
library(scrapeR)
library(rARPACK)
library(Matrix)
library(DT)
library(scales)
#library(FlashGraphR)

source('utils.R')


options(shiny.maxRequestSize = 8000*1024^2)

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

load("openConnectomeGraphs.RData")

invariants <- data.frame(Invariant=c(
                           "Order",
                           "Size",
                           "#Components",
                           "Max Component Size",
                           "Max Degree",
                           "Average Degree",
                           "Min Degree",
                           "Degree Variance",
                           "Diameter",
                           "Girth",
                           "Density",
                           "Clique Number",
                           "Average Path Length",
                           "Centralization Degree",
                           "Reciprocity",
                           "Transitivity",
                           "Degree Assortativity",
                           "Maximal Coreness",
                           "Fast Greedy Communities",
                           "Greedy Domination Number",
                           "Power Law Fit (Exponent)"),
                        Function=c(
                           'vcount',
                           'ecount',
                           'no.clusters',
                           'max.component.size',
                           'get.max.degree',
                           'get.avg.degree',
                           'get.min.degree',
                           'get.var.degree',
                           'diameter',
                           'get.girth',
                           'graph.density',
                           'clique.number',
                           'average.path.length',
                           'get.cent',
                           'get.rec',
                           'transitivity',
                           'get.core',
                           'assortativity.degree',
                           'get.fastgreedy',
                           'dominate.num.greedy',
                           'get.pl'),
                       stringsAsFactors=FALSE)
                           

                         
