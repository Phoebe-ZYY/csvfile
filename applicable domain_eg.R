#applicable domain


library(applicable)
data(qsar_binary)
binary_tr
jacc_sim <- apd_similarity(binary_tr)
jacc_sim

library(ggplot2)

# Plot the empirical cumulative distribution function for the training set
autoplot(jacc_sim)

# Summarize across all training set similarities
mean_sim <- score(jacc_sim, new_data = binary_unk)
mean_sim



#We will use the Ames IA housing data for our example.

library(AmesHousing)
ames <- make_ames()

library(recipes)
library(dplyr)

# Load custom houses from applicable.
data(ames_new, package = "applicable")

ames_cols <- names(ames_new)

training_data <- 
  ames %>% 
  # For consistency, only analyze the data on new properties
  dplyr::select(one_of(ames_cols)) %>% 
  mutate(
    # There is a new neighborhood in ames_new
    Neighborhood = as.character(Neighborhood),
    Neighborhood = factor(Neighborhood, levels = levels(ames_new$Neighborhood))
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
ames_pca <- apd_pca(training_recipe, training_data, threshold = 0.25)
ames_pca



library(ggplot2)
autoplot(ames_pca)

autoplot(ames_pca, matches("PC0[1-5]"))
autoplot(ames_pca, distance) + scale_x_log10()

ames_pca <- apd_pca(training_recipe, training_data)
pca_score <- score(ames_pca, ames_new)
pca_score %>% select(matches("PC00[1-9]"), contains("distance"))

training_scores <- score(ames_pca, training_data)
ggplot(training_scores, aes(x = PC010)) + 
  geom_histogram(col = "white", binwidth = .5) + 
  geom_vline(xintercept = pca_score$PC010, col = "red")


