################################################################################
### k=5, vMF mixture, S^2000 - experiment
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
  message(rep)
  set.seed(rep)
  d <- 2001
  c1 <- rep(0,d); c1[1] <- 1
  c2 <- rep(0,d); c2[2] <- 1
  c3 <- rep(0,d); c3[3] <- 1
  c4 <- rep(0,d); c4[1] <- -1
  c5 <- rep(0,d); c5[2] <- -1
  count <- table(sample(c(1,1,2,2,3,3,4,5),200,replace=TRUE))
  vmf5 <- rvmf(count[1],c1,600)
  vmf5 <- rbind(vmf5, rvmf(count[2],c2,600))
  vmf5 <- rbind(vmf5, rvmf(count[3],c3,600))
  vmf5 <- rbind(vmf5, rvmf(count[4],c4,600))
  vmf5 <- rbind(vmf5, rvmf(count[5],c5,600))
  
  distmat_vmf5 <- matrix(0,ncol=200,nrow=200)
  for(i in 1:199){
    for(j in (i+1):200){
      distmat_vmf5[i,j] <- DIST(vmf5[i,],vmf5[j,])
    }
  }
  distmat_vmf5 <- distmat_vmf5 + t(distmat_vmf5)
  
  res_vmf5 <- mccluster_stability(vmf5, reduced_dim=401, B1=100, delta=0.1, 
                                  train_num=130, tol=0.95, k_max=9, method="skmeans")
  final_result[rep,1] <- which.min(colMeans(res_vmf5)) + 1
  final_result[rep,2] <- (NbClust(vmf5, distance = NULL, diss = distmat_vmf5, max.nc = 10, method = "kmeans", index = "gap")$Best.nc)[1]
  final_result[rep,3] <- (NbClust(vmf5, distance = NULL, diss = distmat_vmf5, max.nc = 10, method = "kmeans", index = "ch")$Best.nc)[1]
  final_result[rep,4] <- (NbClust(vmf5, distance = NULL, diss = distmat_vmf5, max.nc = 10, method = "kmeans", index = "kl")$Best.nc)[1]
  final_result[rep,5] <- (NbClust(vmf5, distance = NULL, diss = distmat_vmf5, max.nc = 10, method = "kmeans", index = "silhouette")$Best.nc)[1]
  final_result[rep,6] <- (NbClust(vmf5, distance = NULL, diss = distmat_vmf5, max.nc = 10, method = "kmeans", index = "hartigan")$Best.nc)[1]
}
final_result