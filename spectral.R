###################################
### Spectral Clust. Function 
###################################

# spectral clustering to k partitions for data on sphere
# X: n x p data matrix
# k: number of clusters
# nn: number of neighbors to use when forming a laplacian matrix
sph_spec_clust <- function(X,k,nn=5){
  N <- nrow(X)
  n.e <- k-1 #how many dimensions to use when applying k-means to the reduced dimension.
  distmat <- matrix(0,ncol=N,nrow=N)
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      distmat[i,j] <- acos(sum(X[i,]*X[j,])) 
    }
  }
  distmat <- distmat + t(distmat)

  knn_mat <- matrix(0,nrow=N,ncol=N)
  for (i in 1:N) {
    neighbor_index <- order(distmat[i,])[2:(nn+1)]
    knn_mat[i,neighbor_index] <- 1 
  }
  knn_mat <- knn_mat + t(knn_mat)
  knn_mat[knn_mat == 2] = 1
  W <- knn_mat 
  
  L <- diag(rowSums(W)) - W
  e <- eigen(L)
  res <- kmeans(e$vectors[,(N-n.e):(N-1)],k)
  return(res)
}
