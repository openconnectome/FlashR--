library(igraph)
library(shinyBS)
library(rARPACK)
library(Matrix)
library(DT)

options(shiny.maxRequestSize = 30*1024^2)

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
                  'Laplacian')

invariants <- c("Degree Distribution",
                "Alpha Centrality",
                "Betweenness",
                "Closeness",
                "Eigenvalue Centrality",
                "Boncich Power Centrality",
                "Authority Score",
                "Hub Score",
                "Articulation Points",
                "Diversity",
                "Average KNN Degree",
                "Page Rank")

