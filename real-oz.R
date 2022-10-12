################################################################################
### Real data - The Oz Books
################################################################################
library(here); here()

source("functions.R")
library(tm)
library(corpus)
library(latex2exp)

# Download the text files, 
# trim header and footer added by the project Guterberg.
# file names denote the author by their first letter, sorted by date of publication
## B=Frank Baum
## T=Ruth Plumly Thompson

filenames <- c(paste0("B",formatC(1:14,width=2,flag=0)),paste0("T",c(15:19,29:33)))
urls <- c("https://www.gutenberg.org/files/55/55-0.txt",             #B01
          "https://www.gutenberg.org/cache/epub/54/pg54.txt",        #B02
          "https://www.gutenberg.org/cache/epub/33361/pg33361.txt",  #B03
          "https://www.gutenberg.org/cache/epub/22566/pg22566.txt",  #B04
          "https://www.gutenberg.org/cache/epub/26624/pg26624.txt",  #B05
          "https://www.gutenberg.org/cache/epub/517/pg517.txt",      #B06
          "https://www.gutenberg.org/cache/epub/955/pg955.txt",      #B07
          "https://www.gutenberg.org/cache/epub/52176/pg52176.txt",  #B08
          "https://www.gutenberg.org/files/51263/51263-0.txt",       #B09
          "https://www.gutenberg.org/cache/epub/25581/pg25581.txt",  #B10
          "https://www.gutenberg.org/cache/epub/959/pg959.txt",      #B11
          "https://www.gutenberg.org/cache/epub/30852/pg30852.txt",  #B12
          "https://www.gutenberg.org/files/50194/50194-0.txt",       #B13
          "https://www.gutenberg.org/cache/epub/961/pg961.txt",      #B14
          "https://www.gutenberg.org/cache/epub/30537/pg30537.txt",  #T15
          "https://www.gutenberg.org/files/53765/53765-0.txt",       #T16
          "https://www.gutenberg.org/cache/epub/58765/pg58765.txt",  #T17
          "https://www.gutenberg.org/files/61681/61681-0.txt",       #T18
          "https://www.gutenberg.org/files/65849/65849-0.txt",       #T19
          "https://www.gutenberg.org/cache/epub/55851/pg55851.txt",  #T29
          "https://www.gutenberg.org/cache/epub/56073/pg56073.txt",  #T30
          "https://www.gutenberg.org/cache/epub/56079/pg56079.txt",  #T31
          "https://www.gutenberg.org/cache/epub/56085/pg56085.txt",  #T32
          "https://www.gutenberg.org/cache/epub/55806/pg55806.txt")  #T33
dir.create(file.path("Oz")) #create a new folder
for(i in 8:length(filenames)){ #trim header and footer of each text and save it locally
  raw <- readLines(urls[i], encoding = "UTF-8")
  start <- grep("^\\*\\*\\* START OF TH", raw) + 1
  if(length(start)==0) start <- grep("^\\*\\*\\*START OF TH", raw) + 1
  stop <- grep("^\\*\\*\\* END OF TH", raw) - 1
  if(length(stop)==0) stop <- grep("^\\*\\*\\*END OF TH", raw) - 1
  lines <- raw[start:stop]
  write(lines, paste0("Oz/",filenames[i],".txt")) 
}

# make DTM (document term matrix) using `tm` package
raw <- VCorpus(DirSource("./Oz", encoding = "UTF-8"))
DTM_raw <- DocumentTermMatrix(raw, control = list(removePunctuation = TRUE, stopwords = TRUE))
sphere_raw <- as.matrix(DTM_raw)
dim(sphere_raw) #24 x 24468
sphere_raw <- sphere_raw / sqrt(rowSums(sphere_raw^2))
mcperturb_accept_rate(sphere_raw, reduced_dim=201, B1=100, delta=0.1, tol=0.95)
set.seed(0)
res_raw <- mccluster_stability(sphere_raw, reduced_dim=201, B1=100, delta=0.1, 
                             train_num=16, k_max=9, tol=0.95, seed=0)

range_all = apply(res_raw[,1:8],2,cumsum)/1:100
par(mfrow=c(1,2))
plot(c(101,101,101,102.5,101,98,101,104),colMeans(res_raw),pch=as.character(2:9),
     xlim=c(1,103), ylim=c(min(range_all),max(range_all)),
     xlab="Number of perturbed samples B", ylab="Instability")
matlines(1:100, apply(res_raw[,1:8],2,cumsum)/1:100, pch=20, col="black",lty=1)
plot(2:9,colMeans(res_raw),type="b",xlab="Number of clusters k",ylab="Instability",pch=20) #select 2

# 6*12 plot



######### how many correctly partitioned? ######################################

set.seed(0)
res_Oz_mat <- matrix(ncol=100,nrow=24)
for(i in 1:100){
  res_Oz_2 <- skmeans(sphere_raw, k=2, control = list(nruns=1))
  res_Oz_mat[,i] <- res_Oz_2$cluster 
}
sapply(1:100, \(i) var(res_Oz_mat[1:14,i]) == 0 && var(res_Oz_mat[15:24,i]) == 0) |> sum()



########## additional experiments for varying (d,epsilon) ######################

for(d in c(200,300,400)){
  for(epsilon in c(0.1,0.2)){
    message(paste0("d=",d,", epsilon=",epsilon))
    print(mcperturb_accept_rate(X=sphere_raw, reduced_dim=d+1, B1=100, delta=epsilon, tol=0.95))
  }
}

res_raw <- list()
i <- 1
for(d in c(200,300,400)){
  for(epsilon in c(0.1,0.2)){
    res_raw[[i]] <- mccluster_stability(sphere_raw, reduced_dim=d+1, B1=100, delta=epsilon, 
                                   train_num=16, k_max=9, tol=0.95, seed=0)
    i <- i+1
  }
}

par(mfrow=c(3,2))

i <- 1
for(d in c(200,300,400)){
  for(epsilon in c(0.1,0.2)){
    plot(2:9,colMeans(res_raw[[i]]),type="b",xlab="Number of clusters k",ylab="Instability",pch=20) 
    title(TeX(paste0("d=",d,", $\\epsilon$=$",epsilon,"$")))
    i <- i+1
  }
} #12*12 plot


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
