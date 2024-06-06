library(Rcpp)

# Define the C++ code with Rcpp
sourceCpp("repid_flex.cpp")

search_split = function(x, y, splits = unique(x)) {
  checkmate::assert_vector(x)
  checkmate::assert_matrix(y)
  
  # Call the C++ function
  objectives = search_split_cpp(x, y, splits)
  
  best = which.min(objectives)
  split.points = splits[best]
  
  best_split_info = find_best_split(x, split.points)
  right = best_split_info$right
  left = best_split_info$left
  
  data.frame(split.points = (left + right) / 2, objective.value = objectives[best])
}

best_split = function(X, y) {
  checkmate::assert_data_frame(X)
  
  # Compute unique quantiles in R and pass to C++ function
  #splits_list = lapply(X, function(x) unique(x))
  
  # Call the C++ function
  res = best_split_cpp(X, y)#, splits_list)
  return(res)
}

# Example
library(data.table)
library(ranger)
library(iml)
library(tidyverse)

set.seed(1)

# Simulate Data
n = 10000
x1 = round(runif(n, -1, 1), 1)
x2 = round(runif(n, -1, 1), 3)
x3 = sample(c(0, 1), size = n, replace = TRUE, prob = c(0.5, 0.5))
x4 = sample(c(0, 1), size = n, replace = TRUE, prob = c(0.7, 0.3))

# noisy vars
x5 = sample(c(0, 1), size = n, replace = TRUE, prob = c(0.5, 0.5))
x6 = rnorm(n, mean = 1, sd = 5)

# target
y = 0.2*x1 - 8*x2 + ifelse(x3 == 0, I(16*x2),0) + ifelse(x1 > 0, I(8*x2),0)
eps = rnorm(n, 0, 0.1*sd(y))
y = y + eps

dat = data.frame(x1, x2, x3, x4, x5, x6, y)
X = dat[, setdiff(colnames(dat), "y")]

# Fit model and compute ICE for x2
mod = ranger(y ~ ., data = dat, num.trees = 100)
pred = function(model, newdata) predict(model, newdata)$predictions
model = Predictor$new(mod, data = X, y = dat$y, predict.function = pred)
effect = FeatureEffect$new(model, method = "ice", grid.size = 50, feature = "x2")

# Center ICE curves
eff = as.data.table(effect$results)
eff = as.data.frame(eff[, .value := (.value - mean(.value)), by = c(".type", ".id")])

# Plot ICE curves: WE WANT TO FIND SUBGROUPS SUCH THAT ICE CURVES ARE HOMOGENOUS
ggplot(eff, aes(x = x2, y = .value)) + 
  geom_line(aes(group = .id))

# Get ICE values and arrange them in a horizontal matrix
Y = spread(eff, x2, .value)
Y = Y[, setdiff(colnames(Y), c(".type", ".id"))]
Y = as.matrix(Y)

str(X) # contains our feature values
str(Y) # contains ICE values for each grid point


microbenchmark::microbenchmark(best_split2(X, Y), best_split(X, Y), times = 50)
