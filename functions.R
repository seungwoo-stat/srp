library(MASS)
library(progress)
library(skmeans)

library(foreach)
library(doParallel)
library(doRNG)

#' Geodesic distance between two points. 
#' 
#' Vectorized operations are not supported.
#' @param a A vector of length 1.
#' @param b A vector of length 1.
#' @return The geodesic distance between \code{a} and \code{b}. 
#' @export
#' @examples
#' DIST(c(1,0,0), c(0,1,0))
DIST <- function(a, b){
  s <- sum(a*b)
  if(s>=1) return(0)
  if(s<=-1) return(pi)
  acos(sum(a*b))
}

#' Returns a logarithmic map
#' 
#' @param p A vector on ambient space of S^d.
#' @param v A vector on ambient space of S^d. 
#' @return Log_p(v) in R^{d+1}. 
#' @export
#' @examples
#' set.seed(0)
#' x <- runif(5)
#' y <- -sqrt((1-x^2)/3)
#' z <- sqrt((1-x^2)/3*2)
#' s1 <- sphere2(x,y,z)
#' LOG(s1[1,],s1[2,])
LOG <- function(p,v){
  d <- DIST(p,v)
  if(d < 1e-10) return(c(0,0,0))
  if(abs(d-pi)<1e-10) stop("log mapping not defined")
  dir <- v - p * cos(d)
  map <- d * dir/sqrt(sum(dir^2))
  return(map)
}

#' Returns an exponential map
#' 
#' @param p A 3-vector on S^d.
#' @param v A 3-vector on T_p(S^d). Thus, p must be orthogonal to v. 
#' @return Exp_p(v). 
#' @export
#' @examples
#' set.seed(0)
#' x <- runif(5)
#' y <- -sqrt((1-x^2)/3)
#' z <- sqrt((1-x^2)/3*2)
#' s1 <- sphere2(x,y,z)
#' EXP(s1[1,],c(0,0,0))
#' EXP(s1[1,],c(-s1[1,2],s1[1,1],0))
EXP <- function(p,v){
  d <- sqrt(sum(v^2))
  if(d < 1e-10) return(p)
  res <- as.vector(p*cos(d) + v/d*sin(d))
  return(res)
}

#' Parallel transport a vector on the tangent space via geodesic
#' 
#' @param start A vector on S^d.
#' @param end A vector on S^d. 
#' @param vec A vector on T_start(S^d), orthogonal to start vector
#' @return Parallel transport \code{vec} from \code{start} to \code{end} via a geodesic. 
#' @export
PARALLEL <- function(start, end, vec){
  if(sum((start+end)^2) < 1e-7) stop("parallel transport to antipodals not supported")
  theta <- DIST(start, end)
  if(sum(theta^2)< 1e-7) return(vec)
  e4 <- LOG(start, end)
  e4 <- e4 / sqrt(sum(e4^2))
  a <- sum(e4*vec) # projection of vec to e4 direction
  stationary <- vec - a*e4
  res <- stationary + a*(e4*cos(theta) - start*sin(theta))
  return(res)
}

#===============================================================================

# cluster stability validation for spherical k-means clustering
# 
# X: data matrix of n x p
# reduced_dim: target dimension of spherical random projection (d in the manuscript)
# B1: number of repetitions (B in the manuscript)
# delta: additive error (epsilon in the manuscript)
# train_num: number of training sets among n data points (m in the manuscript)
# k_max: evaluate number of clusters for 2:k_max
# tol: tolerance, usually set to 0.95
cluster_stability <- function(X, reduced_dim, B1, delta, train_num, k_max, tol=1){
  N <- nrow(X)
  dim <- ncol(X)
  ori_dist <- matrix(0,ncol=N,nrow=N)
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      ori_dist[i,j] <- DIST(X[i,],X[j,])
    }
  }
  ori_dist <- ori_dist + t(ori_dist)
  instability <- matrix(ncol=k_max-1,nrow=B1)
  colnames(instability) <- 2:k_max
  pb <- progress_bar$new(format = " Progress: [:bar] :percent, Estimated completion time: :eta",
                         total = B1)
  
  for(b in 1:B1){
    pb$tick()
    # Random matrix (dim x reduced_dim)
    R1 <- matrix(rnorm(dim*reduced_dim),nrow=dim, ncol=reduced_dim)
    R2 <- matrix(rnorm(dim*reduced_dim),nrow=dim, ncol=reduced_dim)
    red_data1 <- X %*% R1 / sqrt(rowSums((X %*% R1)^2))
    red_data2 <- X %*% R2 / sqrt(rowSums((X %*% R2)^2))
    red_dist1 <- matrix(0,ncol=N,nrow=N)
    red_dist2 <- matrix(0,ncol=N,nrow=N)
    
    while(TRUE){
      red_dist1 <- matrix(0,ncol=N,nrow=N)
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          red_dist1[i,j] <- DIST(red_data1[i,],red_data1[j,])
        }
      }
      red_dist1 <- red_dist1 + t(red_dist1)
      if(mean(red_dist1 <= ori_dist+delta & red_dist1 >= ori_dist-delta)>=tol) break()
      R1 <- matrix(rnorm(dim*reduced_dim),nrow=dim, ncol=reduced_dim)
      red_data1 <- X %*% R1 / sqrt(rowSums((X %*% R1)^2))
    }
    while(TRUE){
      red_dist2 <- matrix(0,ncol=N,nrow=N)
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          red_dist2[i,j] <- DIST(red_data2[i,],red_data2[j,])
        }
      }
      red_dist2 <- red_dist2 + t(red_dist2)
      if(mean(red_dist2 <= ori_dist+delta & red_dist2 >= ori_dist-delta)>=tol) break()
      R2 <- matrix(rnorm(dim*reduced_dim),nrow=dim, ncol=reduced_dim)
      red_data2 <- X %*% R2 / sqrt(rowSums((X %*% R2)^2))
    }
    
    I <- sort(sample(1:N,train_num))  #training index
    J <- setdiff(1:N,I)               #testing index
    for(kc in 2:k_max){
      res1 <- skmeans(red_data1[I,],k=kc,control=list(nruns=5))
      res2 <- skmeans(red_data2[I,],k=kc,control=list(nruns=5))
      centroid1 <- res1$prototypes
      centroid2 <- res2$prototypes
      test1_membership <- vector(length = length(J))
      test2_membership <- vector(length = length(J))
      for(j in 1:length(J)){
        test1_membership[j] <- which.min(acos(colSums(red_data1[J[j],] * t(centroid1))))
        test2_membership[j] <- which.min(acos(colSums(red_data2[J[j],] * t(centroid2))))
      }
      dist_clust <- 0
      for(l1 in 1:(length(J)-1)){
        for(l2 in l1:length(J)){
          dist_clust <- dist_clust+abs((test1_membership[l1]==test1_membership[l2])-
                                         (test2_membership[l1]==test2_membership[l2]))
        }
      }
      dist_clust <- dist_clust/(length(J)*(length(J)-1)/2)
      instability[b,kc-1] <- dist_clust
    }
  }
  return(instability)
}

# cluster stability validation for spherical k-means clustering
# uses multiple cores, on default, uses (maximum number of cores - 2)
mccluster_stability <- function(X, reduced_dim, B1, delta, train_num, k_max, tol=1, seed=123){
  N <- nrow(X)
  dim <- ncol(X)
  ori_dist <- matrix(0,ncol=N,nrow=N)
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      ori_dist[i,j] <- DIST(X[i,],X[j,])
    }
  }
  ori_dist <- ori_dist + t(ori_dist)

  cl <- parallel::makeCluster(detectCores()-2, outfile = "")
  registerDoParallel(cl)
  pb <- txtProgressBar(min = 1, max = B1, style = 3)
  instability <- foreach(b=1:B1, .packages='skmeans', .combine = rbind, .options.RNG=seed) %dorng% {
    setTxtProgressBar(pb, b) 
    instab_temp <- vector(length=k_max-1)
    
    # Random matrix (dim x reduced_dim)
    R1 <- matrix(rnorm(dim*reduced_dim),nrow=dim, ncol=reduced_dim)
    R2 <- matrix(rnorm(dim*reduced_dim),nrow=dim, ncol=reduced_dim)
    red_data1 <- X %*% R1 / sqrt(rowSums((X %*% R1)^2))
    red_data2 <- X %*% R2 / sqrt(rowSums((X %*% R2)^2))
    red_dist1 <- matrix(0,ncol=N,nrow=N)
    red_dist2 <- matrix(0,ncol=N,nrow=N)
    
    while(TRUE){
      red_dist1 <- matrix(0,ncol=N,nrow=N)
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          red_dist1[i,j] <- acos(sum(red_data1[i,]*red_data1[j,]))
        }
      }
      red_dist1 <- red_dist1 + t(red_dist1)
      if(mean(red_dist1 <= ori_dist+delta & red_dist1 >= ori_dist-delta)>=tol) break()
      R1 <- matrix(rnorm(dim*reduced_dim),nrow=dim, ncol=reduced_dim)
      red_data1 <- X %*% R1 / sqrt(rowSums((X %*% R1)^2))
    }
    while(TRUE){
      red_dist2 <- matrix(0,ncol=N,nrow=N)
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          red_dist2[i,j] <- acos(sum(red_data2[i,]*red_data2[j,]))
        }
      }
      red_dist2 <- red_dist2 + t(red_dist2)
      if(mean(red_dist2 <= ori_dist+delta & red_dist2 >= ori_dist-delta)>=tol) break()
      R2 <- matrix(rnorm(dim*reduced_dim),nrow=dim, ncol=reduced_dim)
      red_data2 <- X %*% R2 / sqrt(rowSums((X %*% R2)^2))
    }
    
    I <- sort(sample(1:N,train_num))  #training index
    J <- setdiff(1:N,I)               #testing index
    for(kc in 2:k_max){
      res1 <- skmeans(red_data1[I,],k=kc,control=list(nruns=5))
      res2 <- skmeans(red_data2[I,],k=kc,control=list(nruns=5))
      centroid1 <- res1$prototypes
      centroid2 <- res2$prototypes
      test1_membership <- vector(length = length(J))
      test2_membership <- vector(length = length(J))
      for(j in 1:length(J)){
        test1_membership[j] <- which.min(acos(colSums(red_data1[J[j],] * t(centroid1))))
        test2_membership[j] <- which.min(acos(colSums(red_data2[J[j],] * t(centroid2))))
      }
      dist_clust <- 0
      for(l1 in 1:(length(J)-1)){
        for(l2 in l1:length(J)){
          dist_clust <- dist_clust+abs((test1_membership[l1]==test1_membership[l2])-
                                         (test2_membership[l1]==test2_membership[l2]))
        }
      }
      dist_clust <- dist_clust/(length(J)*(length(J)-1)/2)
      instab_temp[kc-1] <- dist_clust
    }
    instab_temp
  }
  close(pb)
  stopCluster(cl)
  colnames(instability) <- 2:k_max
  return(instability)
}

############################

# before conducting the *cluster_stability function, run this function to determine the reduced_dim
# X: data matrix of n x p
# reduced_dim: target dimension of spherical random projection (d in the manuscript)
# B1: number of repetitions (B in the manuscript)
# delta: additive error (epsilon in the manuscript)
# tol: tolerance, usually set to 0.95
perturb_accept_rate <- function(X, reduced_dim, B1, delta, tol=1){
  N <- nrow(X)
  dim <- ncol(X)
  ori_dist <- matrix(0,ncol=N,nrow=N)
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      ori_dist[i,j] <- DIST(X[i,],X[j,])
    }
  }
  ori_dist <- ori_dist + t(ori_dist)

  acceptance_rate <- vector(length=B1)
  pb <- progress_bar$new(format = " Progress: [:bar] :percent, Estimated completion time: :eta",
                         total = B1)
  for(b in 1:B1){
    pb$tick()
    # Random matrix (dim x reduced_dim)
    R1 <- matrix(rnorm(dim*reduced_dim),nrow=dim, ncol=reduced_dim)
    red_data1 <- X %*% R1 / sqrt(rowSums((X %*% R1)^2))

    # acceptance rate
    red_dist1 <- matrix(0,ncol=N,nrow=N)
    for(i in 1:(N-1)){
      for(j in (i+1):N){
        red_dist1[i,j] <- DIST(red_data1[i,],red_data1[j,])
      }
    }
    red_dist1 <- red_dist1 + t(red_dist1)
    acceptance_rate[b] <- (mean(red_dist1 <= ori_dist+delta & red_dist1 >= ori_dist-delta)>=tol)
  }
  return(mean(acceptance_rate))
}

# perturb_accept_rate function with multiple cores
# on default, uses (maximum number of cores - 2)
mcperturb_accept_rate <- function(X, reduced_dim, B1, delta, tol=1, seed=123){
  N <- nrow(X)
  dim <- ncol(X)
  ori_dist <- matrix(0,ncol=N,nrow=N)
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      ori_dist[i,j] <- DIST(X[i,],X[j,])
    }
  }
  ori_dist <- ori_dist + t(ori_dist)
  
  acceptance_rate <- vector(length=B1)
  cl <- parallel::makeCluster(detectCores()-2, outfile = "")
  registerDoParallel(cl)
  pb <- txtProgressBar(min = 1, max = B1, style = 3)
  acceptance_rate <- foreach (b=1:B1, .combine=c, .verbose = F, .options.RNG=seed) %dorng% {
    # pb$tick()
    setTxtProgressBar(pb, b) 
    # Random matrix (dim x reduced_dim)
    R1 <- matrix(rnorm(dim*reduced_dim),nrow=dim, ncol=reduced_dim)
    red_data1 <- X %*% R1 / sqrt(rowSums((X %*% R1)^2))
    
    # acceptance rate
    red_dist1 <- matrix(0,ncol=N,nrow=N)
    for(i in 1:(N-1)){
      for(j in (i+1):N){
        red_dist1[i,j] <- acos(sum(red_data1[i,]*red_data1[j,]))
      }
    }
    red_dist1 <- red_dist1 + t(red_dist1)
    mean(red_dist1 <= ori_dist+delta & red_dist1 >= ori_dist-delta)>=tol
  }
  close(pb)
  stopCluster(cl)
  return(mean(acceptance_rate))
}

############################

# cluster stability validation for spectral clustering
# uses multiple cores, on default, uses (maximum number of cores - 2)
mccluster_stability_spec <- function(X, reduced_dim, B1, delta, train_num, k_max, tol=1, seed=123){
  N <- nrow(X)
  dim <- ncol(X)
  ori_dist <- matrix(0,ncol=N,nrow=N)
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      ori_dist[i,j] <- DIST(X[i,],X[j,])
    }
  }
  ori_dist <- ori_dist + t(ori_dist)
  
  cl <- parallel::makeCluster(detectCores()-2, outfile = "")
  registerDoParallel(cl)
  pb <- txtProgressBar(min = 1, max = B1, style = 3)
  instability <- foreach(b=1:B1, .combine = rbind, .options.RNG=seed) %dorng% {
    source("spectral.R")
    setTxtProgressBar(pb, b) 
    instab_temp <- vector(length=k_max-1)
    
    # Random matrix (dim x reduced_dim)
    R1 <- matrix(rnorm(dim*reduced_dim),nrow=dim, ncol=reduced_dim)
    R2 <- matrix(rnorm(dim*reduced_dim),nrow=dim, ncol=reduced_dim)
    red_data1 <- X %*% R1 / sqrt(rowSums((X %*% R1)^2))
    red_data2 <- X %*% R2 / sqrt(rowSums((X %*% R2)^2))
    red_dist1 <- matrix(0,ncol=N,nrow=N)
    red_dist2 <- matrix(0,ncol=N,nrow=N)
    
    while(TRUE){
      red_dist1 <- matrix(0,ncol=N,nrow=N)
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          red_dist1[i,j] <- acos(sum(red_data1[i,]*red_data1[j,]))
        }
      }
      red_dist1 <- red_dist1 + t(red_dist1)
      if(mean(red_dist1 <= ori_dist+delta & red_dist1 >= ori_dist-delta)>=tol) break()
      R1 <- matrix(rnorm(dim*reduced_dim),nrow=dim, ncol=reduced_dim)
      red_data1 <- X %*% R1 / sqrt(rowSums((X %*% R1)^2))
    }
    while(TRUE){
      red_dist2 <- matrix(0,ncol=N,nrow=N)
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          red_dist2[i,j] <- acos(sum(red_data2[i,]*red_data2[j,]))
        }
      }
      red_dist2 <- red_dist2 + t(red_dist2)
      if(mean(red_dist2 <= ori_dist+delta & red_dist2 >= ori_dist-delta)>=tol) break()
      R2 <- matrix(rnorm(dim*reduced_dim),nrow=dim, ncol=reduced_dim)
      red_data2 <- X %*% R2 / sqrt(rowSums((X %*% R2)^2))
    }
    
    I <- sort(sample(1:N,train_num))  #training index
    J <- setdiff(1:N,I)               #testing index
    nn.choose <- 10                   #nn to choose the membership
    for(kc in 2:k_max){
      res1 <- sph_spec_clust(red_data1[I,],k=kc)
      res2 <- sph_spec_clust(red_data2[I,],k=kc)
      train1_membership <- res1$cluster
      train2_membership <- res2$cluster
      test1_membership <- vector(length = length(J))
      test2_membership <- vector(length = length(J))
      for(j in 1:length(J)){
        nn1_memberships <- train1_membership[which(rank(acos(colSums(red_data1[J[j],] * t(red_data1[I,])))) <= nn.choose)]
        nn2_memberships <- train2_membership[which(rank(acos(colSums(red_data2[J[j],] * t(red_data2[I,])))) <= nn.choose)]
        ux <- unique(nn1_memberships)
        test1_membership[j] <- ux[which.max(tabulate(match(nn1_memberships, ux)))]
        ux <- unique(nn2_memberships)
        test2_membership[j] <- ux[which.max(tabulate(match(nn2_memberships, ux)))]
      }
      dist_clust <- 0
      for(l1 in 1:(length(J)-1)){
        for(l2 in l1:length(J)){
          dist_clust <- dist_clust+abs((test1_membership[l1]==test1_membership[l2])-
                                         (test2_membership[l1]==test2_membership[l2]))
        }
      }
      dist_clust <- dist_clust/(length(J)*(length(J)-1)/2)
      instab_temp[kc-1] <- dist_clust
    }
    instab_temp
  }
  close(pb)
  stopCluster(cl)
  colnames(instability) <- 2:k_max
  return(instability)
}