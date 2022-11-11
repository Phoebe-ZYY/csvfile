## This script is the first part for reproducing results for 
##  1NN tangent canformal method. The result is a file with p-values.
start0 <- Sys.time() # timing start
## 1. Load data :
usps.train <- as.matrix(read.table(file = normalizePath("USPS//USPS.train"), skip = 3))
usps.test <- as.matrix(read.table(file = normalizePath("USPS//USPS.test"), skip = 3))
n.train <- nrow(usps.train)
n.test <- nrow(usps.test)
usps <- rbind(usps.train, usps.test)

## normalise objects
x <- usps[, 1 : 256]
x.mean <- rowMeans(x)
x <- x - x.mean
deviation <- rowSums(x^2) / ncol(x)
x <- x / sqrt(deviation)
usps[, 1 : 256] <- x

## 2. Shuffle data 
seed <- 0
set.seed(0)
n <- n.train + n.test
indices <- sample(n)
permutted <- usps[indices,]
x.train <- permutted[1 : n.train, -257]
y.train <- permutted[1 : n.train, 257]
x.test <- permutted[(n.train + 1) : n, -257]

## 3. Calculate tangent distances:
## The package TangentDistance provides an interface for 
## Daniel Keyser's implementation of tangent distance.
## The package can be downloaded from 
##  http://www.cs.rhul.ac.uk/~valentina/R-packages/TangentDistance_1.0.tar.gz
## To install the package uncomment the following two lines correcting the path to the dowloaded archive:
# path <- "TangentDistance//Rpackage//TangentDistance_1.0.tar.gz"
# install.packages( normalizePath(path), repos = NULL, type = "source")

library(TangentDistance)
start <- Sys.time() # timing for tangent distances
tangent.dist.matrix.train <- tangentDistMatrix(xtrain = x.train, 
                                                height = 16, width = 16)
end <- Sys.time() # timing for tangent distances end
time <- difftime(end, start, units = "mins") 
cat("train distances elapsed time", time, " minutes\n")

start <- Sys.time() # timing for tangent distances
tangent.dist.matrix.test <- tangentDistMatrix(xtrain = x.train,
                                               xtest = x.test, 
                                               height = 16, width = 16) 
end <- Sys.time() # timing for tangent distances end
time <- difftime(end, start, units = "mins") 
cat("test distances elapsed time", time, " minutes\n")

## 4. Calculate conformity scores using 1-NN method
labels <- 0 : 9
n.labels <- length(labels)
# First find nearest neighbours for training examples to optimise the calculation.
cat("sort distances for training objects\n")
train.my.dist <- rep(Inf, n.train)
train.other.dist <- rep(Inf, n.train)
train.conf.scores <- rep(0, n.train)
for(i in 1 : n.train ) {
  if(i%%100==0) { #just to show that the algorithm is calculating
    cat("|")
    if(i%%1000==0)  cat("\n",i) 
  }
  #distances to other objects
  if((i > 1) && (i < n.train)) {
    col.indices <- 1 : (i - 1)
    # distances to previous objects are in the row (i - 1) 
    tmp1<- tangent.dist.matrix.train[(col.indices - 1) * n.train - 
                         col.indices * (col.indices - 1) / 2 + (i - col.indices)]
    # distances to objects after the current object are in the column i, rows ((i + 1) : n.train) 
    row.indices <- (i + 1) : n.train 
    tmp2 <- tangent.dist.matrix.train[(i - 1) * n.train - 
                          i * (i - 1) / 2 + (row.indices - i)]
    dist.vector <- c(tmp1, tmp2)
    
  } else if(i == 1) {
    # for the first train object
    dist.vector <- tangent.dist.matrix.train[ 1 : (n.train - 1)]
  } else if(i == n.train) {
    # for the last train object
    col.indices <- 1 : (n.train - 1)
    dist.vector <- tangent.dist.matrix.train[(col.indices - 1) * n.train - 
                                 col.indices * (col.indices - 1) / 2 + (i - col.indices)]
  }
  my.dist <- dist.vector[y.train[-i] == y.train[i]]
  other.dist <- dist.vector[y.train[-i] != y.train[i]]
  train.my.dist[i] <- min(my.dist) 
  train.other.dist[i] <- min(other.dist)
  train.conf.scores[i] <- train.other.dist[i] / train.my.dist[i]
}
cat("\n") 

cat("calculate conformity scores\n")
less.conform.test <- matrix(nrow = n.test, ncol = n.labels)
eq.conform.test <- matrix(nrow = n.test, ncol = n.labels)
for(i in 1 : n.test) {
  if(i%%100==0) { #just to show that the algorithm is calculating
    cat("|")
    if(i%%1000==0)  cat("\n",i) 
  }
  #distances 
  dist.vector <- tangent.dist.matrix.test[i, ]
  
  # try labels
  for(j in 1 : n.labels) {
    conf.scores <- train.conf.scores
    
    cur.label <- labels[j]
    my.dist <- min(dist.vector[y.train == cur.label])
    other.dist <- min(dist.vector[y.train != cur.label])
    test.conf.score <- other.dist / my.dist 

    #train conformity scores
    for(i.train in 1 : n.train) {
      if(y.train[i.train] == cur.label) {
       if(dist.vector[i.train] < train.my.dist[i.train])
        conf.scores[i.train] <- train.other.dist[i.train] / dist.vector[i.train] 
      } 
    }
    my.idx <- y.train == cur.label
    n.my <- sum(my.idx) + 1 # including the test example
    my.scores <- conf.scores[my.idx]
    
    lessN <- sum(my.scores < test.conf.score)
    eqN <- sum(my.scores == test.conf.score)  + 1
    
    less.conform.test[i, j] <- lessN / n.my
    eq.conform.test[i, j] <- eqN / n.my
    
  }
}
cat("\n")


## 5. Calculate p-values
set.seed(0) # for generating \tau for p-values
pvalues <- matrix(nrow = n.test, ncol = n.labels)
colnames(pvalues) <- labels
taus <- runif(n.test)
for(i in 1 : n.test) {
  for(j in 1 : n.labels){
    pvalues[i, j] <- less.conform.test[i, j] +
                       taus[i] * eq.conform.test[i, j]
  }
}
  
file.result <- normalizePath(
  paste("usps_pvalues_1NN_tangent_seed=", seed,".txt", sep = ""), 
  mustWork = FALSE)
write.table(x = pvalues, file = file.result, 
            col.names = TRUE, row.names = FALSE)
end0 <- Sys.time() # timing end
time <- difftime(end0, start0, units = "mins") 
cat("elapsed time", time, "minutes\n")
