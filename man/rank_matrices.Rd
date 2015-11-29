\name{rank_matrices}
\alias{rank_matrices}
\title{Rank Matrices}
\usage{
  ranked <- rank_matrices(graphs, normalize)
  
  Inputs:
    graphs[N rois x N rois x k scans]: a matrix containing 1 NxN graph for each subject in the dataset
    normalize[boolean]: an option to normalize the ranked graphs and fix their values between 0 and 1
    Outputs:
    ranked[N rois x N rois x k scans]: a matrix where each of the k scans are ranked internally using the built in rank function    
}
\description{
The Rank Matrices function ranks each scan of a set of graphs internally based on their value. Ranked graphs have been shown to afford better reliability than unranked graphs, and provide an alternative path for statistical assessment that will serve value for future data explorations.
}