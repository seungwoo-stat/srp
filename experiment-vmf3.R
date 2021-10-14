################################################################################
### k=3, vMF mixture, S^1000 - experiment
################################################################################

library(here); here()

source("functions.R")

library(Directional) #rvmf
library(latex2exp)
library(RColorBrewer)
library(skmeans) #spk-means
library(NbClust)

final_result <- matrix(nrow=50,ncol=6)
colnames(final_result) <- c("proposed","gap","ch","kl","silhouette","hartigan")

for(rep in 1:50){
  set.seed(rep)
  d <- 1001
  c1 <- rep(0,d); c1[1] <- 1
  c2 <- rep(0,d); c2[2] <- 1
  c3 <- rep(0,d); c3[3] <- 1
  count <- table(sample(1:3,150,replace=TRUE))
  vmf3 <- rvmf(count[1],c1,300)
  vmf3 <- rbind(vmf3, rvmf(count[2],c2,300))
  vmf3 <- rbind(vmf3, rvmf(count[3],c3,300))
  
  distmat_vmf3 <- matrix(0,ncol=150,nrow=150)
  for(i in 1:149){
    for(j in (i+1):150){
      distmat_vmf3[i,j] <- DIST(vmf3[i,],vmf3[j,])
    }
  }
  distmat_vmf3 <- distmat_vmf3 + t(distmat_vmf3)

  mcperturb_accept_rate(vmf3, reduced_dim=201, B1=100, delta=0.2, tol=0.99)
  res_vmf3 <- mccluster_stability(vmf3, reduced_dim=201, B1=100, delta=0.2, 
                                  train_num=100, tol=0.99, k_max=9, method="skmeans")
  final_result[rep,1] <- which.min(colMeans(res_vmf3)) + 1
  
  final_result[rep,2] <- (NbClust(vmf3, distance = NULL, diss = distmat_vmf3, max.nc = 10, method = "kmeans", index = "gap")$Best.nc)[1]
  final_result[rep,3] <- (NbClust(vmf3, distance = NULL, diss = distmat_vmf3, max.nc = 10, method = "kmeans", index = "ch")$Best.nc)[1]
  final_result[rep,4] <- (NbClust(vmf3, distance = NULL, diss = distmat_vmf3, max.nc = 10, method = "kmeans", index = "kl")$Best.nc)[1]
  final_result[rep,5] <- (NbClust(vmf3, distance = NULL, diss = distmat_vmf3, max.nc = 10, method = "kmeans", index = "silhouette")$Best.nc)[1]
  final_result[rep,6] <- (NbClust(vmf3, distance = NULL, diss = distmat_vmf3, max.nc = 10, method = "kmeans", index = "hartigan")$Best.nc)[1]
}
final_result
