
library(reshape2)
library(stringr)
library(ggplot2)
library(tidyverse)
library(Rmalschains)
library(iml)
library(ranger)
library(kmlShape)
library(dtw)
library(tidyr)

#------------------------------------------------------------------------------------------------------------
# helper functions

generate_split_candidates = function(xval, use.quantiles = TRUE) {
  if (use.quantiles) { # to speedup we use only quantile values as possible split points
    q = unique(quantile(xval, seq(0.01, 0.99, by = 0.01), type = 1))
  } else {
    q = sort(unique(xval))
  }
  # use mid-point of two subsequent split candidates
  q = q + c(diff(q)/2, 0)
  
  if (max(q) == max(xval))
    q = q[-length(q)]
  
  return(q)
}


# Frechet distance measure (sum of squares)
SS_fre = function(y, x, requires.x = FALSE) { # slow
  # using only y-axis of curves is enough as x-axis is always the same for all curves
  require(kmlShape)
  center = colMeans(y)
  grid.x = as.numeric(names(center))
  pdp.y = unname(center)
  dist = apply(y, 1, function(ice) distFrechet(grid.x, pdp.y, grid.x, ice, FrechetSumOrMax = "sum"))
  sum(dist)
}


# NEW:
# Filter function to use only ice curves within grid points to find best split point
filter = function(feat, X, Y, sub_number){
  values = unique(X[sub_number,feat])
  dt <- data.table(value = values)
  dt[, Range := cut(values, breaks = c(min(X[,feat])-0.0001,as.numeric(names(Y)), max(X[,feat])+0.0001), labels = c(names(Y), names(Y)[length(names(Y))]))]
  return(which(names(Y) %in% unique(dt$Range)))
}



#------------------------------------------------------------------------------------------------------------
# adjuat perform.split and optimization function:

# adjust perform_split: filter y.list - use only ice curves for which data points of feature combination is 
# available, otherwise throw them out --> no extrapolation -> calculate objective only on these curves and find
# with that best split point
perform_split = function(split.points, feat, x, xval, y, min.node.size, objective) {

  
  split.points = sort.int(split.points)
  # assign intervalnr. according to split points
  node.number = findInterval(x = xval, split.points, rightmost.closed = TRUE) + 1
  
  # compute size of each childnode
  node.size = tabulate(node.number)
  # if minimum node size is violated, return Inf
  if (min(node.size) < min.node.size)
    return(Inf)
  # split y according to node.number
  y.list = split(y, node.number)
  # use only ice curves that are available for the combination of the two features -> no extrapolation
  y.list.filtered = sapply(seq_along(unique(node.number)), FUN = function(i) {
    sub_number = which(node.number==i)
    ind = filter(feat, x, y, sub_number)
    y.list[[i]] = y.list[[i]][,ind]
  }, USE.NAMES = FALSE)
  # x.list only needed if this is used in the objective
  requires.x = formals(objective)[["requires.x"]]
  if (isTRUE(requires.x))
    x.list = split(xval, node.number) else
      x.list = NULL
  
  res = vapply(seq_along(y.list.filtered), FUN = function(i) {
    objective(y = y.list.filtered[[i]], x = x.list[[i]])
  }, FUN.VALUE = NA_real_, USE.NAMES = FALSE)
  list(sum(res), y.list.filtered)
}


# optimization function: test with binary split - return also ice curves for plotting
find_best_binary_split = function(feat, x, xval, y, n.splits = 1, min.node.size = 1, 
                                  objective, ...) {
  assert_choice(n.splits, choices = 1)
  
  # use different split candidates to perform split
  q = generate_split_candidates(xval, use.quantiles = TRUE)
  splits = vapply(q, FUN = function(i) {
    perform_split(i, feat = feat, x = x, xval = xval, y = y, min.node.size = min.node.size, 
                  objective = objective)[[1]]
  }, FUN.VALUE = NA_real_, USE.NAMES = FALSE)
  # select the split point yielding the minimal objective
  best = which.min(splits)
  # return also ice curves of best split point for plotting
  ice = perform_split(q[best], feat = feat, x = x, xval = xval, y = y, min.node.size = min.node.size, 
                objective = objective)[[2]]
  
  return(list(split.points = q[best], objective.value = splits[best], ice.curves = ice))
}


#------------------------------------------------------------------------------------------------------------
# Simulation Study

# Simulate Data
set.seed(1234)

z <- runif(500,0,1)
d = ifelse(z>0.5,1,0)
eps = rnorm(500, 0, 0.15)

# noisy vars
x5 = sample(c(0, 1), size = 500, replace = TRUE, prob = c(0.5, 0.5))
x6 = rnorm(500, mean = 1, sd = 5)


y = z^0.5 * (1-d) + (2^(1.5*z)-1) * d + eps
dat = data.frame(z,d,x5,x6,y)
X = dat[, setdiff(colnames(dat), "y")]
#write.csv2(dat, "data/ExtrapolationExample1.csv")


# Fit model and compute ICE for z
mod = ranger(y ~ ., data = dat, num.trees = 500)
pred = predict.function = function(model, newdata) predict(model, newdata)$predictions
model = Predictor$new(mod, data = X, y = dat$y, predict.function = pred)
effect = FeatureEffects$new(model, method = "ice", grid.size = 20, features = "z")
# Get ICE values and arrange them in a horizontal matrix
Y = spread(effect$results$z, .borders, .value)
Y = Y[, setdiff(colnames(Y), c(".type", ".id", ".feature"))]


# Plot ICE curves
ggplot(effect$results$z, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id))

#------------------------------------------------------------------------------------------------------------

# Results with adjusted optimization function of binary splits:

# get results of best binary splits and ice curves of best split point
res = find_best_binary_split(feat = "z", x = X, xval = X[,"d"], y = Y, objective = SS_fre)
plot.data.new = res$ice.curves

# prepare data for plotting
plot.prep = function(i, x){
  x$.id = 1:nrow(x)
  x = gather(x, .borders, .value, colnames(x)[1]:colnames(x)[ncol(x)-1], factor_key=TRUE)
  x$.split = i
  return(x)
}

data.new = lapply(seq_along(plot.data.new), FUN = function(i) {
  plot.prep(i, plot.data.new[[i]])
})
data.new = rbindlist(data.new)


# Plot splitted ICE curves
ggplot(data.new, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id)) + facet_grid(~ .split)
