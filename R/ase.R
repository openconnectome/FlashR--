ase <- function(A, dims){
  if(nrow(A) >= 400){
    require(irlba)
    A.svd <- irlba(A, nu = dims, nv = dims)
    A.svd.values <- A.svd$d[1:dims]
    A.svd.vectors <- A.svd$v[,1:dims]
    if(dim == 1)
      A.coords <- sqrt(A.svd.values) * A.svd.vectors
    else
      A.coords <- A.svd.vectors %*% diag(sqrt(A.svd.values))
  } else{
    A.svd <- svd(A)
    if(dims == 1)
      A.coords <- A.svd$v[,1] * sqrt(A.svd$d[1])
    else
      A.coords <- A.svd$v[,1:dims] %*% diag(sqrt(A.svd$d[1:dims]))
  }

  return(A.coords)
}
