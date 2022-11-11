## This script is the second part for reproducing results for 
##  1NN tangent canformal method. This part uses p-values calculated in the first part.
start0 <- Sys.time() # timing start
## 1. Load data :
usps.train <- as.matrix(read.table(file = normalizePath("USPS//USPS.train"), skip = 3))
usps.test <- as.matrix(read.table(file = normalizePath("USPS//USPS.test"), skip = 3))
n.train <- nrow(usps.train)
n.test <- nrow(usps.test)
usps <- rbind(usps.train, usps.test)

## 2. Shuffle data
seed <- 0
set.seed(seed)
n <- n.train + n.test
indices <- sample(n)
permutted <- usps[indices,]
y.test <- permutted[(n.train + 1) : n, 257]
## 3. Load p-values
file.result <- normalizePath(
  paste("usps_pvalues_1NN_tangent_seed=", seed,".txt", sep = ""), 
  mustWork = TRUE)
pvalues <- read.table(file = file.result, header = TRUE)
labels <- 0 : 9
n.labels <- length(labels)
## 4. Convert p-values into probabilities
  prob <- matrix(nrow = n.test, ncol = n.labels)
  colnames(prob) <- labels 
  library("fdrtool") # for implementation of Grenander's estimator
  for(j in 1 : n.labels) {
    gr.res <- grenander(ecdf(c(0, pvalues[, j], 1)))
    f <- gr.res$f.knots
    x <- gr.res$x.knots
    n <- length(f)
    D <- f[n]

    for(i in 1 : n.test){
      p <- pvalues[i, j]
      if((p < x[n]) && (p > x[1]) ) 
        f.j <- f[ which(x > p)[1] - 1 ] 
      else if (p >= x[n]) # it means that p = 1
        f.j <-  f[ n ] 
      else # it means that p = 0
        f.j <-  f[ 1 ]
      prob[i, j] <- D / f.j
    }
      
  }
  prob.sum <- rowSums( prob )
  
  prob.norm <- matrix(nrow = n.test, ncol = n.labels)
  for(i in 1 : n.test)
   prob.norm[i,] <- prob[i,] / prob.sum[i]
  colnames(prob.norm) <- labels  

## 5. Calculate losses
# log loss  
log.loss <- 0
for (i in 1 : n.test) 
  log.loss <- log.loss - log(prob.norm[i, as.character(y.test[i])]) 
log.loss <- log.loss / n.test
# square loss
sq.loss <- 0
for (i in 1 : n.test) {
  v1 <- rep(0, n.labels)
  v1[y.test[i] + 1] <- 1
  v2 <- prob.norm[i, ]
  sq.loss <- sq.loss + sum((v1 - v2)^2) 
}  
sq.loss <- sqrt(sq.loss / (n.labels * n.test))
 
## Results
res.string <- paste("1NN tangent conformal probabilities: log loss =", log.loss,
                     ", square loss =", sq.loss, "\n")
cat(res.string)
end0 <- Sys.time() # timing end
time <- difftime(end0, start0, units = "mins") 
cat("elapsed time", time, "minutes\n")
