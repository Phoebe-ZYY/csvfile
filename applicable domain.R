#applicable domain

application <- application_imputation_combi %>%
  ungroup() %>%
  select(-name_ai) %>%
  select(-name_crop_earthstat) %>%
  select(-name_country) %>%
  select(-casrn_chemical) %>%
  distinct() %>%
  filter(!is.na(landevap)) %>%
  filter(!is.na(HDI))
  
application_n <- application_n_imputation_combi %>%
  ungroup() %>%
  select(-name_ai) %>%
  select(-name_crop_earthstat) %>%
  select(-name_country) %>%
  select(-casrn_chemical) %>%
  distinct()%>%
  filter(!is.na(landevap)) %>%
  filter(!is.na(HDI))

library(applicable)
library(ggplot2)
library(recipes)
library(dplyr)


app_cols <- names(application_n)

levels(factor(application$layer_climate))

training_data <- 
  application %>% 
  # For consistency, only analyze the data on new properties
  dplyr::select(one_of(app_cols)) %>% 
  mutate(
    # There is a new climate layer EF in application_n
    layer_climate = as.character(layer_climate),
    layer_climate = factor(layer_climate, levels = c(levels(factor(application$layer_climate)), "EF")),
    # There is 2 new moa group H30 I30 in application_n
    name_moagroup = as.character(name_moagroup),
    name_moagroup = factor(name_moagroup, levels = c(levels(factor(application$name_moagroup)), "H30", "I30"))
  )


training_recipe <-
  recipe( ~ ., data = training_data) %>%
  step_dummy(all_nominal()) %>% 
  # Remove variables that have the same value for every data point.
  step_zv(all_predictors()) %>% 
  # Transform variables to be distributed as Gaussian-like as possible.
  step_YeoJohnson(all_numeric()) %>%
  # Normalize numeric data to have a mean of zero and
  # standard deviation of one.
  step_normalize(all_numeric())

ames_pca <- apd_pca(training_recipe, training_data)
ames_pca

#Since no threshold was provided, the function computed the number of 
#principal components that accounted for at most 95% of the total variance.

#For illustration, setting threshold = 0.25 or 25%, 
#we now need only 10 principal components:
#ames_pca <- apd_pca(training_recipe, training_data, threshold = 0.25)
#ames_pca


autoplot(ames_pca)

autoplot(ames_pca, matches("PC0[1-5]"))
autoplot(ames_pca, distance) + scale_x_log10()

ames_pca <- apd_pca(training_recipe, training_data)
pca_score <- score(ames_pca, application_n)
pca_score %>% select(matches("PC00[1-9]"), contains("distance"))

mean(pca_score$distance) #12.4
mean(pca_score$distance_pctl)

training_scores <- score(ames_pca, training_data)

training_scores$flag <- "training"
pca_score$flag <- "pca"
score <- rbind(training_scores, pca_score)

ggplot(score, aes(x = PC001, fill = flag)) +
  geom_histogram(binwidth = .5, alpha=.5, position="identity")
  
# `ames_pca$pcs` is the output of `prcomp()`
comp_one <- ames_pca$pcs$rotation[, 2]
comp_one[order(abs(comp_one), decreasing = TRUE)] %>% head(5)

summary(ames_pca$pcs)
ames_pca.df <- as.data.frame(ames_pca$pcs$rotation)

training_data$flag <- "training"
application_n$flag <- "pca"
application_data <- rbind(training_data, application_n)

ggplot(application_data, aes(x = name_indication, fill = flag)) + 
  geom_histogram(binwidth = .5, alpha=.5, position="identity", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(pca_score, aes(x = distance)) +
  geom_histogram(binwidth = .5, alpha=.5, position="identity")

n_distinct(c(application$name_cropgroup, application$name_chemclass1, 
             application$name_moagroup, application$name_indication,
             application$layer_climate))
