ase <- function(A, dim){
  if(nrow(A) >= 400){
    require(irlba)
    A.svd <- irlba(A, nu = dim, nv = dim)
    A.svd.values <- A.svd$d[1:dim]
    A.svd.vectors <- A.svd$v[,1:dim]
    if(dim == 1)
      A.coords <- sqrt(A.svd.values) * A.svd.vectors
    else
      A.coords <- A.svd.vectors %*% diag(sqrt(A.svd.values))
  } else{
    A.svd <- svd(A)
    if(dim == 1)
      A.coords <- A.svd$v[,1] * sqrt(A.svd$d[1])
    else
      A.coords <- A.svd$v[,1:dim] %*% diag(sqrt(A.svd$d[1:dim]))
  }
  
  return(A.coords)
}