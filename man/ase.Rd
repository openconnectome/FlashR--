\name{ase}
\alias{ase}
\title{Adjacency Spectral Embedding}
\usage{
A_embed <- ase(A, dim)
}
\description{
The eigendecomposition of an adjacency matrix provides a way to embed a graph  as points in finite dimensional Euclidean space.  This embedding allows the full arsenal of statistical and machine learning methodology  for multivariate Euclidean data to be deployed for graph inference. Adjacency Spectral Embedding performs this embedding of matrix A into dim dimensions.
}
\examples{
N <- 200
A <- matrix(data = runif(N*N), ncol = N, nrow = N)
dims <- 15
A_embed <- ase(A, dims)
}
