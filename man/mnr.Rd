\name{mnr}
\alias{mnr}
\title{Mean Normalized Rank}
\usage{
me_rank <- mnr(rdf, remove_outliers, thresh, output)

Inputs:
  rdf[(N scans*k scans/subject) x 1]: the reliability distribution function; see rdf documentation for details
  remove_outliers[boolean]: indicates whether outliers below the threshold are to be removed from the calculation; defaults to TRUE
  thresh[double]: the threshold for removing outliers
  output[boolean]: indicates whether to print a summary of the removed outliers is to be printed
Outputs:
  me_rank[double]: a scalar indicating the mean normalized rank
}
\description{
The mean normalized rank, or mnr, of a set of scans gives an indication as to the reliability of each subject. A high mnr corresponds to lower intra subject distances and high reliability, whereas a low mnr corresponds to higher intra subject distances and correspondingly low reliablity. If intra subject relationships are indistinguishable from inter subject relationships, the graphs are said to go to chance and will have an mnr approaching .5. 
}
