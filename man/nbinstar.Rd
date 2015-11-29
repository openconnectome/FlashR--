\name{nbinstar}
\alias{nbinstar}
\title{Nbin*}
\usage{
nbinpack <- nbinstar(graphs, ids, scans, N, spacing, lim)

Inputs:
  graphs[N rois x N rois x k scans double]: the set of graphs for each scan
  ids[N scans string]: a vector of strings, with the subject id of each scan
  scans[integer]: the number of scans for each subject; default is 2
  N[integer]: the number of bins to test; default is 100
  spacing['linear' or 'log']: the spacing to place between the bins; default is 'linear'
  lim[integer]: the maximum number of bins to test; default is 0
Outputs:
  king[double]: the highest mnr value during the nbin thresholding procedure
  kingbins[integer]: the number of bins at the optimal mnr
  kingrdf[(N scans*kscans/subject)x1 double]: a vector indicating the reliability distsribution function of the dataset at the highest mnr
  mnrs[N bins tested x 1 double]: the list of mnrs associated with each number of bins
  set[N bins tested x 1 integer]: the list of bins tested
}
\description{
The nbinstar function thresholds graphs based on a fixed number of bins. For example, if we select 5 bins with linear spacing, the values under .2 will receive a value of .2, the values between .2 and .4 will receive a value of .4, and so on. We then compute the mnr at each of these number of bins, and look for the optimal number of bins to maximize the mnr. nbins = N^2 can be said to be a fully ranked situation, where all of the values for each graph, each subject are ranked from lowest to greatest.
}
