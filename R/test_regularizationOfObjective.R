source("R/perform_split.R")
source("R/split_optimizer.R")
source("R/split_parent_node.R")


library(reshape2)
library(stringr)
library(ggplot2)
library(tidyverse)
library(Rmalschains)
library(iml)
library(ranger)
library(kmlShape)
library(dtw)



#---------------------------------------------------------------------------------------------------------------
# Frechet distance measure (sum of squares)
SS_fre = function(y, x, requires.x = TRUE) { # slow
  # using only y-axis of curves is enough as x-axis is always the same for all curves
  require(kmlShape)
  center = colMeans(y)
  grid.x = as.numeric(names(center))
  pdp.y = unname(center)
  dist = apply(y, 1, function(ice) distFrechet(grid.x, pdp.y, grid.x, ice, FrechetSumOrMax = "sum"))
  sum(dist)
}


#---------------------------------------------------------------------------------------------------------------
# Simulation Example

# Simulate Data
set.seed(1234)
n = 500
x1 = runif(n, -1, 1)
x2 = runif(n, -1, 1)
x3 = sample(c(0, 1), size = n, replace = TRUE, prob = c(0.5, 0.5))
x4 = sample(c(0, 1), size = n, replace = TRUE, prob = c(0.7, 0.3))
eps = rnorm(n, 0, 1)

# noisy vars
x5 = sample(c(0, 1), size = n, replace = TRUE, prob = c(0.5, 0.5))
x6 = rnorm(n, mean = 1, sd = 5)

# new y_hat with conitnuous interactions
y = 0.2*x1 - 8*x2^2 + 5*cos(x2*5)*x6 + ifelse(x3 == 0, I(8*x2),0) + eps
dat = data.frame(x1, x2, x3, x4, x5, x6, y)
X = dat[, setdiff(colnames(dat), "y")]


# Fit model 
mod = ranger(y ~ ., data = dat, num.trees = 1000)
pred = predict.function = function(model, newdata) predict(model, newdata)$predictions
model = Predictor$new(mod, data = X, y = dat$y, predict.function = pred)



#---------------------------------------------------------------------------------------------------------------

# 1. Suggestion to adjust objective:

# Idea: if there is an interaction, then subgroup PDPs differ a lot from each other and from the overall PDP
# --> regularization of objective by calculating distance of clusterPDPs to overallPDP
# Steps:
# 1. calculate within_SS (as before with SS_fre)
# 2. calcuale between_SS by calculating frechetdistance between cluster PDPs and overall PDP
# 3. calculate new objective: sum(within_SS) + lambda * (within_SS / (between_SS/#clusters))

# since large values of between_SS mean that there is more interaction and vice versa, we need to divide by this
# value since we have a mnimization problem. To regularize also for the number of clusters, we divide also by #clusters
# difficulty: how to choose lambda?

# perform.split function was adjusted. Additional parameters are dist.between (choose between this method "fre" or
# "sd" - see below) and lambda 




# x2 interacts with x6 and x3, results with this method for lambda = 0.25 for multiway splits

# compute ICE for x2
effect = FeatureEffects$new(model, method = "ice", grid.size = 20, features = "x2")

# Plot ICE curves
ggplot(effect$results$x2, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id))

# Get ICE values and arrange them in a horizontal matrix
Y = spread(effect$results$x2, .borders, .value)
Y = Y[, setdiff(colnames(Y), c(".type", ".id", ".feature"))]

# find best multiway split for parend node 
SS_fre(Y, X) #(38766)
sp_multi = lapply(1:5, function(i) {
  split_parent_node(Y = Y, X = X, objective = SS_fre, dist.between = "fre", lambda = 0.25,
                    n.splits = i, optimizer = find_best_multiway_split2, min.node.size = 10)
})
# results show, that x6 reduces original SS by ~8000, also x3 reduces it by more than 3000.
# Due to the regularization, objective vaules can be higher than the overall SS of 38766
# Number of splits are a bit regularized e.g. as can be seen for x2, but not too much



# x1 does not interact with any other feature, 
# results with this method for lambda = 0.25 for multiway splits

# compute ICE for x1
effect = FeatureEffects$new(model, method = "ice", grid.size = 20, features = "x1")

# Plot ICE curves
ggplot(effect$results$x1, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id))

# Get ICE values and arrange them in a horizontal matrix
Y = spread(effect$results$x1, .borders, .value)
Y = Y[, setdiff(colnames(Y), c(".type", ".id", ".feature"))]

# find best multiway split for parend node 
SS_fre(Y, X) #(61577.32)
sp_multi = lapply(1:3, function(i) {
  split_parent_node(Y = Y, X = X, objective = SS_fre, dist.between = "fre", lambda = 0.25,
                    n.splits = i, optimizer = find_best_multiway_split2, min.node.size = 10)
})

sp_multi
# [[1]]
# feature objective.value split.points best.split
# 1:      x1        61687.40   -0.7243768      FALSE
# 2:      x2        55803.31   -0.5138996       TRUE
# 3:      x3        63352.57          0.5      FALSE
# 4:      x4        61953.10          0.5      FALSE
# 5:      x5        62453.05          0.5      FALSE
# 6:      x6        56072.02     11.41767      FALSE

# although we regularize, x2 still reduces the objectve funtion by 5600 and almost the same for x6 for a binary split
# in fact the value did not much increase to the objective without regularization which was 55713 for x2.
# Also if we increase lambda to 0.5, it increases only to 55893. --> it seams the cluster pdps here to actually
# differ to overall PDP

# Visualize results e.g. for binary split
node_index_frechet = generate_node_index(Y, X, result = sp_multi[[1]])
plot.data = effect$results$x1
plot.data$.split = node_index_frechet$class[plot.data$.id]

ggplot(plot.data, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id)) + facet_grid(~ .split)



#---------------------------------------------------------------------------------------------------------------

# 2. Suggestion to adjust objective:

# Idea: if there is an interaction, then VARIATION of subgroup PDPs differ a lot from each other and from the 
# VARIATION of the overall PDP
# --> regularization of objective by calculating deviation of sd of clusterPDPs to sd of overallPDP
# Steps:
# 1. calculate within_SS (as before with SS_fre)
# 2. calcuale SD of cluster PDPs and overall PDP
# 3. calculate SD_Dev = abs(weighted.mean(sd_clusterPDP) - sd_overallPDP)
# 3. calculate new objective: sum(within_SS) + lambda * (within_SS / SD_Dev)

# since large values of SD_DEV mean that there is more interaction and vice versa, we need to divide by this
# value since we have a mnimization problem. 
# We use weighted mean to adjust for high variations within clusters with few observations
# difficulty: how to choose lambda?


# x2 interacts with x6 and x3, results with this method for lambda = 0.015 for multiway splits

# compute ICE for x2
effect = FeatureEffects$new(model, method = "ice", grid.size = 20, features = "x2")

# Plot ICE curves
ggplot(effect$results$x2, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id))

# Get ICE values and arrange them in a horizontal matrix
Y = spread(effect$results$x2, .borders, .value)
Y = Y[, setdiff(colnames(Y), c(".type", ".id", ".feature"))]

# find best multiway split for parent node 
SS_fre(Y, X) #(38766)
sp_multi = lapply(1:5, function(i) {
  split_parent_node(Y = Y, X = X, objective = SS_fre, dist.between = "sd", lambda = 0.015,
                    n.splits = i, optimizer = find_best_multiway_split2, min.node.size = 10)
})
sp_multi
# results show, that x6 reduces original SS by ~6000, also x3 reduces it a bit (not too much, eventually smaller value for lambda more suitable)
# Due to the regularization, objective vaules can be higher than the overall SS of 38766 --> larger effect than "fre"
# Number of splits are not regularized 



# x1 does not interact with any other feature, 
# results with this method for lambda = 0.015 for multiway splits

# compute ICE for x1
effect = FeatureEffects$new(model, method = "ice", grid.size = 20, features = "x1")

# Plot ICE curves
ggplot(effect$results$x1, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id))

# Get ICE values and arrange them in a horizontal matrix
Y = spread(effect$results$x1, .borders, .value)
Y = Y[, setdiff(colnames(Y), c(".type", ".id", ".feature"))]

# find best multiway split for parend node 
SS_fre(Y, X) #(61577.32)
sp_multi = lapply(1:3, function(i) {
  split_parent_node(Y = Y, X = X, objective = SS_fre, dist.between = "sd", lambda = 0.015,
                    n.splits = i, optimizer = find_best_multiway_split2, min.node.size = 10)
})

sp_multi
# [[1]]
# feature objective.value split.points best.split
# 1:      x1       264704.31    0.7523979      FALSE
# 2:      x2        64200.61   -0.4997017       TRUE
# 3:      x3        67826.75          0.5      FALSE
# 4:      x4       115120.81          0.5      FALSE
# 5:      x5        79736.05          0.5      FALSE
# 6:      x6        85411.18    -5.818304      FALSE
# 
# [[2]]
# feature objective.value          split.points best.split
# 1:      x1       108718.83     0.763805,0.825875      FALSE
# 2:      x2        59671.05 -0.4962238, 0.6640457       TRUE
# 3:      x3        67826.75                   0.5      FALSE
# 4:      x4       115120.81                   0.5      FALSE
# 5:      x5        79736.05                   0.5      FALSE
# 6:      x6        68562.42    4.709172,10.516012      FALSE

# none of the features reduces the original objective value of 61577 in a binary split (but values decrease fast
# for more splits e.g. x2 --> how to regularize on that??)
# also for x5 no ojective value of any feature is reduced in the binary split. --> also here the problem of 
# fast decreasing values with more splits 




