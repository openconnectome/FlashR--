\name{rdf}
\alias{rdf}
\title{Reliability Distribution Function}
\usage{
r_dist <- rdf(dist, ids, scans)

Inputs:
  dist [N subjects x N subjects]: a distance matrix containing the graph distances between subjects
  ids [N subjectsx1]: a vector containing the string ids of each subject
  scans [k]: a scalar indicating the maximum number of scans per subject in the dataset
  
Outputs:
  r_dist[(N_subjects*k scans) x 1]: a rank comparison vector between all intra subject scans
}
\description{
The reliability distribution function provides a method of analyzing a distance matrix for a
set of subjects, and analyzing how close the intra subject relationship is. This function allows statistical analysis of the reliability of a dataset.
}
\examples{
  # see documentation for compute_distance for details about the graphs parameter
  dist <- compute_distance(graphs)
  r_dist <- rdf(dist, ids, scans)
}