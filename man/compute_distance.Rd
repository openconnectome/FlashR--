\name{compute_distance}
\alias{compute_distance}
\title{Distance Computation}
\usage{
  dist <- compute_distance(graphs, norm)
  
  Inputs:
    graphs[N rois x N rois x k scans]: a matrix containing 1 NxN graph for each subject in the dataset
    norm (defaults to 'F'): the normalization method for acquiring a scalar from the difference between subjects
  
  Outputs:
    dist[k scans x k scans]: a matrix containing the distance between each subject
}
\description{
  The distance computation allows analyzing the distance between a set of graphs. Each subject
  is compared pairwise, first by computing the distance at each roi by matrix subtraction, and then
  acquiring a scalar from this NxN result by normalizing. The distance matrix allows for analysis
  of the distances between scans in a dataset.
}