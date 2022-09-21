#####################################################
### Bulls eye data set (k=2) - experiment
#####################################################

library(here); here()

source("functions.R")
source("spectral.R")

library(NbClust)

final_result <- matrix(nrow=50,ncol=6)
colnames(final_result) <- c("proposed","gap","ch","kl","silhouette","hartigan")

for(rep in 1:50){
  message(rep)
  set.seed(rep)
  e1 <- c(1,rep(0,1000))
  bull <- matrix(ncol=1001, nrow=400)
  
  for(i in 1:200){
    v1 <- rnorm(3)
    v1 <- pi/8*v1/sqrt(sum(v1^2))
    v2 <- rnorm(3)
    v2 <- pi/4*v1/sqrt(sum(v1^2))
    ep1 <- rnorm(1000)
    ep1 <- pi/16*ep1/sqrt(sum(ep1^2))
    ep2 <- rnorm(1000)
    ep2 <- pi/16*ep2/sqrt(sum(ep2^2))
    # near data
    t1 <- EXP(e1,c(0,v1,rep(0,997)))
    bull[i,] <- EXP(t1,PARALLEL(e1,t1,c(0,ep1)))
    # far data
    t2 <- EXP(e1,c(0,v2,rep(0,997)))
    bull[200+i,] <- EXP(t2,PARALLEL(e1,t2,c(0,ep2)))
  }
  
  N <- 400
  distmat <- matrix(0,ncol=N,nrow=N)
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      distmat[i,j] <- acos(sum(bull[i,]*bull[j,]))
    }
  }
  distmat <- distmat + t(distmat)
  
  res_spec <- mccluster_stability_spec(X=bull, reduced_dim = 151, B1=100, delta=0.1, train_num = 266, k_max=9, tol=0.95)
  
  final_result[rep,1] <- which.min(colMeans(res_spec)) + 1
  final_result[rep,2] <- (NbClust(bull, distance = NULL, diss = distmat, max.nc = 9, method = "kmeans", index = "gap")$Best.nc)[1]
  final_result[rep,3] <- (NbClust(bull, distance = NULL, diss = distmat, max.nc = 9, method = "kmeans", index = "ch")$Best.nc)[1]
  final_result[rep,4] <- (NbClust(bull, distance = NULL, diss = distmat, max.nc = 9, method = "kmeans", index = "kl")$Best.nc)[1]
  final_result[rep,5] <- (NbClust(bull, distance = NULL, diss = distmat, max.nc = 9, method = "kmeans", index = "silhouette")$Best.nc)[1]
  final_result[rep,6] <- (NbClust(bull, distance = NULL, diss = distmat, max.nc = 9, method = "kmeans", index = "hartigan")$Best.nc)[1]
}

final_result

## plot ########################################################################

N <- 400
distmat <- matrix(0,ncol=N,nrow=N)
for(i in 1:(N-1)){
  for(j in (i+1):N){
    distmat[i,j] <- acos(sum(bull[i,]*bull[j,]))
  }
}
distmat <- distmat + t(distmat)

res_1 <- mccluster_stability(X=bull, reduced_dim = 151, B1=100, delta=0.1, train_num = 266, k_max=9, tol=0.95, seed=0)
res_2 <- mccluster_stability_spec(X=bull, reduced_dim = 151, B1=100, delta=0.1, train_num = 266, k_max=9, tol=0.95, seed=0)

par(mfrow=c(1,2))
plot(isoMDS(distmat[1:400,1:400],k=2)$points,pch=ifelse(res_spec$cluster==1,"-","+"),
     xlab="MDS 1", ylab="MDS 2",main="Multidimensional scaling plot")
plot(2:9,colMeans(res_2),type="b",ylim=range(colMeans(res_1)),
     xlab="Number of clusters k",ylab="Instability",pch=20,
     main="Instability plot for two clustering algorithms") 
lines(2:9, colMeans(res_1),lty=2)
points(2:9, colMeans(res_1),pch=20)
legend("topright",lty=1:2,legend=c("spectral","spherical k-means"),bty="n")


