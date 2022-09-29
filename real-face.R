################################################################################
### Real data - The JAFFE Dataset
### * We assume the data set is placed in the `jaffe` folder.
### * We do not contain the data set in this repository.
### * The data set must be obtained via https://zenodo.org/record/3451524#.YWfL8i_kEUE , 
### * with permission from the authors of the data set.
################################################################################
library(here); here()

source("functions.R")
library(tiff)
library(pixmap)
library(latex2exp)

file.name <- list.files("jaffe")
file.name <- file.name[endsWith(file.name,"tiff")]
face <- matrix(0, ncol = 128*128, nrow = 213)
for(i in 1:length(file.name)){
  temp <- rowsum(pixmapGrey(readTIFF(paste0("jaffe/",file.name[i])))@grey, rep(1:128, each=2))/2
  face[i, ] <- rowsum(t(temp), rep(1:128, each=2))/2
}
rownames(face) <- paste0(substr(file.name, 1, 2), "-",substr(file.name, 4,6))
sel_index <- (substr(file.name,1,2)=="NM") | (substr(file.name,1,2)=="KL") | 
  (substr(file.name,1,2)=="KR") | (substr(file.name,1,2)=="UY") | (substr(file.name,1,2)=="YM")
exp_index <- (substr(file.name,4,5)=="HA") | (substr(file.name,4,5)=="NE") | 
  (substr(file.name,4,5)=="SU")
face <- face[sel_index & exp_index,]
name_face <- substr(file.name,1,2)[sel_index & exp_index]
dim(face) # 44 x 16384

file.name[sel_index & exp_index]

draw_index <- c(1,10,18,27,36,4,12,21,30,39,7,15,24,33,42)
par(mfrow=c(3,6))
par(mar=c(1,1,1,1))
for(d in 1:15){
  if(d %% 5 == 1){
    par(mar=c(0,0,0,0))
    plot(0,0,pch="",axes=F)
    text(.75,0,c("HAP","NEU","SUR")[(d %/% 5) +1])
    par(mar=c(1,1,1,1))
  }
  image((matrix(face[draw_index[d],],ncol=128,nrow=128)),col=grey(seq(0, 1, length = 256)),asp=1,ylim=c(1,0),axes=F)
  if(d %/% 6 == 0) title(c("KL","KR","NM","UY","YM")[d %% 6])
}
par(mfrow=c(1,1))
par(mar=c(5.1, 4.1, 4.1, 2.1))

dist_face <- matrix(0,ncol=nrow(face),nrow=nrow(face))
for(i in 1:(nrow(face)-1)){
  for(j in (i+1):nrow(face)){
    dist_face[i,j] = cor(face[i,],face[j,])
  }
}
dist_face <- dist_face + t(dist_face)
isoMDS(dist_face,k=2)$points |> plot(pch=as.numeric(as.factor(name_face)))

sphere_raw <- face / sqrt(rowSums(face^2))

mcperturb_accept_rate(X=sphere_raw, reduced_dim=31, B1=100, delta=0.1, tol=0.95)

set.seed(0)
res_raw <- mccluster_stability(sphere_raw, reduced_dim=201, B1=100, delta=0.1, 
                             train_num=31, k_max=9, tol=0.95, seed=0)

range_all = apply(res_raw[,1:8],2,cumsum)/1:100

par(mfrow=c(1,2))
plot(c(101,101,101,104,101,104,101,104),colMeans(res_raw),pch=as.character(2:9),
     xlim=c(1,103), ylim=c(min(range_all),max(range_all)),
     xlab="Number of perturbed samples B", ylab="Instability")
matlines(1:100, apply(res_raw[,1:8],2,cumsum)/1:100, pch=20, col="black",lty=1)
plot(2:9,colMeans(res_raw),type="b",xlab="Number of clusters k",ylab="Instability",pch=20) #select 5

# 6*12 plot



######### how many correctly partitioned? ######################################

set.seed(0)
res_face_mat <- matrix(ncol=100,nrow=44)
for(i in 1:100){
  res_face_2 <- skmeans(sphere_raw, k=5, control = list(nruns=200))
  res_face_mat[,i] <- res_face_2$cluster 
}

sapply(1:100, \(i) var(res_face_mat[1:9,i]) == 0 && 
         var(res_face_mat[10:17,i]) == 0 &&
         var(res_face_mat[18:26,i]) == 0 &&
         var(res_face_mat[27:35,i]) == 0 &&
         var(res_face_mat[36:44,i]) == 0) |> sum()



########## additional experiments for varying (d,epsilon) ######################

for(d in c(30,50,70)){
  for(epsilon in c(0.1,0.2)){
    message(paste0("d=",d,", epsilon=",epsilon))
    print(mcperturb_accept_rate(X=sphere_raw, reduced_dim=d+1, B1=100, delta=epsilon, tol=0.95))
  }
}

res_raw <- list()
i <- 1
for(d in c(30,50,70)){
  for(epsilon in c(0.1,0.2)){
    res_raw[[i]] <- mccluster_stability(sphere_raw, reduced_dim=d+1, B1=100, delta=epsilon, 
                                        train_num=31, k_max=9, tol=0.95, seed=0)
    i <- i+1
  }
}

par(mfrow=c(3,2))

i <- 1
for(d in c(30,50,70)){
  for(epsilon in c(0.1,0.2)){
    plot(2:9,colMeans(res_raw[[i]]),type="b",xlab="Number of clusters k",ylab="Instability",pch=20) 
    title(TeX(paste0("d=",d,", $\\epsilon$=$",epsilon,"$")))
    i <- i+1
  }
}

# 12* 12 plot

# X <- sphere_raw
# N <- nrow(X)
# dim <- ncol(X)
# ori_dist <- matrix(0,ncol=N,nrow=N)
# for(i in 1:(N-1)){
#   for(j in (i+1):N){
#     ori_dist[i,j] <- DIST(X[i,],X[j,])
#   }
# }
# ori_dist <- ori_dist + t(ori_dist)
# hist(ori_dist[ori_dist != 0],breaks=50,freq=FALSE,xlab="Pairwise distances")
# summary(ori_dist[ori_dist!= 0])