
#application method
## 1. Load data :
result_full_app_p <- read.csv("C:/me/data gaps/20220615result_full_app_p.csv" )
#result_full_app_perc_in <- result_full_app_p %>%
  #filter(!is.na(aerialapplicationp)) %>%
  #filter(!is.na(boomsprayerp)) %>%
  #filter(!is.na(handoperatedsprayerp)) %>%
  #filter(!is.na(tunnelsprayerp))

#scenarios for other 3 application methods without tunnelsprayeryes = 0 (only tunnelsprayeryes has output 0)
result_full_app_p_3 <- result_full_app_p %>%
  filter(tunnelsprayeryes == 0) %>%
  #mutate(boomsprayersudo = 0.01)
  mutate(boomsprayersudo = case_when(boomsprayeryes < 0.5 ~ 0.01,
                                     boomsprayeryes >= 0.5 ~ 0.99))
#pvalues_3 <- result_full_app_p_3[,c(12,16,20)]
pvalues_3 <- result_full_app_p_3[,c(12,37,20)]

#scenarios for all application methods (all yes p-values not equal to zero)
result_full_app_p_4 <- result_full_app_p %>%
  filter(tunnelsprayeryes != 0)
pvalues_4 <- result_full_app_p_4[,c(12,16,20,32)]

#for 3 application methods
labels <- 0:2
n.labels <- length(labels)
n.test <- nrow(result_full_app_p_3)

#for 4 application methods
labels <- 0:3
n.labels <- length(labels)
n.test <- nrow(result_full_app_p_4)

## 2. Shuffle data
seed <- 0
set.seed(seed)
#n <- n.train + n.test
#indices <- sample(n)
#permutted <- usps[indices,]
#y.test <- permutted[(n.train + 1) : n, 257]
## 3. Load p-values

## 4. Convert p-values into probabilities
#should run twice for both with 3 application methods and with 4 application methods
#change pvalues_3 or pvalues_4
prob <- matrix(nrow = n.test, ncol = n.labels)
colnames(prob) <- labels 
library("fdrtool") # for implementation of Grenander's estimator
model <- for(j in 1 : n.labels) {
  gr.res <- grenander(ecdf(c(0, pvalues_3[, j], 1)))
  f <- gr.res$f.knots
  x <- gr.res$x.knots
  n <- length(f)
  D <- f[n]
  
  for(i in 1 : n.test){
    p <- pvalues_3[i, j]
    if((p < x[n]) && (p > x[1]) ) 
      f.j <- f[ which(x > p)[1] - 1 ] 
    else if (p >= x[n]) # it means that p = 1
      f.j <-  f[ n ] 
    else # it means that p = 0
      f.j <-  f[ 1 ]
    prob[i, j] <- D / f.j
  }
  
}
prob.sum <- rowSums(prob)

prob.norm <- matrix(nrow = n.test, ncol = n.labels)
for(i in 1 : n.test){
  prob.norm[i,] <- prob[i,] / prob.sum[i]
  }
colnames(prob.norm) <- labels 
prob.df_3 <- as.data.frame(prob)
prob.norm_3 <- as.data.frame(prob.norm)
cpp_prob_3 <- cbind(pvalues_3, prob, prob.norm)
cpp_prob_3_perc <- cpp_prob_3[,c(1:3,7:9)]
cpp_prob_3_perc$tunnelsprayeryes <- 0
cpp_prob_3_perc$`3` <- 0
cpp_prob_3_perc <- cpp_prob_3_perc %>%
  select(aerialapplicationyes, boomsprayeryes, handoperatedsprayeryes, tunnelsprayeryes, `0`, `1`, `2`, `3`)

cpp_prob_3_perc <- cpp_prob_3_perc %>%
  mutate(aerialperc = `0`, boomperc = `1`, handperc = `2`) %>%
  select(-c(`0`, `1`, `2`))

cpp_prob_3_perc %>%
  select(aerialapplicationyes, aerialperc) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = log(value), fill = name)) +
  geom_histogram(position = "identity", alpha = 0.5)

cpp_prob_3_perc %>%
  select(aerialapplicationyes, aerialperc) %>% 
  ggplot(aes(x = aerialapplicationyes, y = aerialperc)) +
  geom_point()

cpp_prob_3_perc %>%
  select(boomsprayersudo, boomperc) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = log(value), fill = name)) +
  geom_histogram(position = "identity", alpha = 0.5)

cpp_prob_3_perc %>%
  select(boomsprayersudo, boomperc) %>% 
  ggplot(aes(x = boomsprayersudo, y = boomperc)) +
  geom_point()

cpp_prob_3_perc %>%
  select(handoperatedsprayeryes, handperc) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = log(value), fill = name)) +
  geom_histogram(position = "identity", alpha = 0.5)

cpp_prob_3_perc %>%
  select(handoperatedsprayeryes, handperc) %>% 
  ggplot(aes(x = handoperatedsprayeryes, y = handperc)) +
  geom_point()



prob.df_4 <- as.data.frame(prob)
prob.norm_4 <- as.data.frame(prob.norm)
cpp_prob_4 <- cbind(pvalues_4, prob, prob.norm)
cpp_prob_4_perc <- cpp_prob_4[,c(1:4, 9:12)]

cpp_prob_3_perc_full <- cbind(result_full_app_p_3, cpp_prob_3_perc)
cpp_prob_4_perc_full <- cbind(result_full_app_p_4, cpp_prob_4_perc)

cpp_prob_4_perc_full1 <- cpp_prob_4_perc_full1 %>%
  select(-c(cpp_prob_4_perc_full1$X.1))


#cpp_prob_full_perc <- rbind(cpp_prob_3_perc_full, cpp_prob_4_perc_full1)
write.csv(cpp_prob_full_perc, "20220621cpp_prob_app_perc.csv")
write.csv(cpp_prob_4_perc_full, "20220706cpp_prob_4_perc_full.csv")
write.csv(cpp_prob_3_perc_full, "20220708cpp_prob_3_perc_full.csv")
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

#bbch cropstage
result_full_bbch <- read.csv("C:/me/data gaps/20220608result_full_bbch.csv" )
result_full_bbch_in <- na.omit(result_full_bbch)

pvalues <- result_full_bbch_in[,c(12,16,20,32,36,40,44,48,52)]
labels <- 0:8
n.labels <- length(labels)
n.test <- nrow(result_full_bbch_in)

## 4. Convert p-values into probabilities
prob <- matrix(nrow = n.test, ncol = n.labels)
colnames(prob) <- labels 
library("fdrtool") # for implementation of Grenander's estimator
model <- for(j in 1 : n.labels) {
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
for(i in 1 : n.test){
  prob.norm[i,] <- prob[i,] / prob.sum[i]
}
colnames(prob.norm) <- labels 
prob.df <- as.data.frame(prob)
prob.norm <- as.data.frame(prob.norm)
cpp_prob_bbch <- cbind(pvalues, prob, prob.norm)
write.csv(cpp_prob_bbch, "20220608cpp_prob_bbch.csv")



#sensitivity analysis
x <- fast99(model = model, factors = c("aerialapplicationyes", 
                                       "boomsprayeryes", 
                                       "handoperatedsprayeryes", 
                                       "tunnelsprayeryes"), n = 75696,
            q = "qunif", q.arg = list(min = 0.0001, max = 1))

tell(x, y = prob)
plot(x)
print(x)
