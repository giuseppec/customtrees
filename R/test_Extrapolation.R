
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
SS_fre = function(y, x, requires.x = TRUE) { # slow
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
  if(length(unique(X[,feat])) > 2){
    dt <- data.table(value = values)
    dt[, Range := cut(values, breaks = c(min(X[,feat])-0.0001,as.numeric(names(Y)), max(X[,feat])+0.0001), labels = c(names(Y), names(Y)[length(names(Y))]))]
    return(which(names(Y) %in% unique(dt$Range)))
    }
  else if(length(unique(X[,feat])) == length(values)){
    #dt = data.table(value = as.numeric(names(Y)))
    return(c(1:length(names(Y))))
  }
  else if(length(values) == 1 & length(unique(X[,feat]))==2){
    if(values < mean(unique(X[,feat]))){
      #dt = data.table(value = as.numeric(names(Y))[1:round(ncol(Y)/2,0)])
      return(c(1:round(ncol(Y)/2,0)))
    }
    else if(values > mean(unique(X[,feat]))){
      #dt = data.table(value = as.numeric(names(Y))[round(ncol(Y)/2,0):ncol(Y)])
      return(c(round(ncol(Y)/2,0):ncol(Y)))
    }
  }
  
  
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
  y.list.filtered = lapply(seq_along(unique(node.number)), FUN = function(i) {
    sub_number = which(node.number==i)
    ind = filter(feat, x, y, sub_number)
    y.list[[i]][,ind, drop = FALSE]
  })
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
  
  return(list(split.points = q[best], objective.value = splits[best]#, ice.curves = ice
              ))
}



# Optimization with MA-LS Chains: finds best split points and returns value of objective function
find_best_multiway_split2 = function(feat, x, xval, y, n.splits = 1, min.node.size = 1, 
                                     objective, control = NULL, ...) {
  unique.x = length(unique(xval))
  # max. number of splits to be performed must be unique.x-1
  if (n.splits >= unique.x)
    n.splits = unique.x - 1
  # if only one split is needed, we use exhaustive search
  if (n.splits == 1)
    return(find_best_binary_split(feat, x, xval, y, n.splits = n.splits,
                                  min.node.size = min.node.size, objective = objective))
  
  # generate initial population
  xval.candidates = generate_split_candidates(xval, use.quantiles = FALSE)
  pop = replicate(n.splits, sample(xval.candidates, size = min(n.splits*10, unique.x)))
  pop = unique(pop)
  
  # replace maximum number of evaluations depending of number of unique x values
  xval.ncomb = prod(rep(unique.x, n.splits))
  if (is.null(control))
    control = Rmalschains::malschains.control(istep = 100, ls = "sw")
  control$istep = min(control$istep, xval.ncomb)
  
  # possible lower and upper values
  lower = rep(min(xval.candidates), n.splits)
  upper = rep(max(xval.candidates), n.splits)
  
  # optimization function
  .perform_split = function(par)
    perform_split(split.points = par, feat = feat, x = x, xval = xval, y = y, min.node.size = min.node.size,
                  objective = objective)[[1]]
  best = Rmalschains::malschains(.perform_split, lower = lower, upper = upper,
                                 initialpop = pop, verbosity = 0, control = control, ...) #   control = malschains.control(istep = 300, ls = "sw"),
  
  return(list(split.points = sort(best$sol), objective.value = best$fitness))
}



split_parent_node = function(Y, X, feat, n.splits = 1, min.node.size = 1, optimizer, 
                             objective, ...) {
  assert_data_frame(X)
  #assert_choice(target.col, choices = colnames(data))
  assert_integerish(n.splits)
  assert_integerish(min.node.size)
  assert_function(objective, args = c("y", "x", "requires.x"))
  assert_function(optimizer, args = c("xval", "y"))
  
  # find best split points per feature - TODO: adjusted that not splitted after ICE variable - option to choose?
  opt.feature = lapply(X[,-which(colnames(X)==feat)], function(xval) {
    optimizer(feat = feat, x = X, xval = xval, y = Y, n.splits = n.splits, min.node.size = min.node.size, 
              objective = objective, ...)
  })
  
  result = rbindlist(lapply(opt.feature, as.data.frame), idcol = "feature")
  result = result[, .(split.points = list(split.points)), by = c("feature", "objective.value"), with = TRUE]
  result$best.split = result$objective.value == min(result$objective.value)
  #result = result[, best.split := objective.value == min(objective.value)]
  return(result)
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
model = Predictor$new(mod, data = X, y = y, predict.function = pred)
effect = FeatureEffects$new(model, method = "ice", grid.size = 20, features = "x5")
# Get ICE values and arrange them in a horizontal matrix
Y = spread(effect$results$z, .borders, .value)
Y = Y[, setdiff(colnames(Y), c(".type", ".id", ".feature"))]
# centered ICE curves
for(i in 1:nrow(Y)){
  Y[i,] = as.numeric(unname(Y[i,])) - mean(as.numeric(unname(Y[i,])))
}


# Plot ICE curves
ggplot(effect$results$z, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id))

#------------------------------------------------------------------------------------------------------------

# Results with adjusted optimization function of binary splits:

# get results of best binary splits and ice curves of best split point
res = find_best_binary_split(feat = "x2", x = X, xval = X[,"x1"], y = Y, objective = SS_fre)
res2 = split_parent_node(Y = Y, X = X, feat = "x3", optimizer = find_best_multiway_split2, n.splits = 1, objective = SS_fre, min.node.size = 10)
sp(X, Y, "x4", 3, 10, SS_fre, 0.1, 0.1)

sp = function(x, y, feat, n.splits, min.node.size, objective, improve.first.split, improve.n.splits){
  res = list()
  for(i in 1:n.splits) {
    res[[i]] = split_parent_node(Y = y, X = x, feat = feat, n.splits = i, min.node.size = min.node.size, 
                     optimizer = find_best_multiway_split2, objective = objective)
    if(i == 1){
      SS_total = objective(y = y, x = x)
      improve = (SS_total-res[[i]][which(res[[i]]$best.split==TRUE), "objective.value"])/SS_total
      if(improve < improve.first.split) return(NULL)
    }
    if(i > 1){
      obj1 = res[[i-1]][which(res[[i-1]]$best.split==TRUE), "objective.value"]
      obj2 = res[[i]][which(res[[i-1]]$best.split==TRUE), "objective.value"]
      improve = (obj1-obj2)/obj1
      if(improve < improve.n.splits) return(list(res[(i-1)], SS_total))
    }
  }
  return(list(res[[i]], SS_total))
}

library(dplyr)
find_potential_interactions <- function(x, model,  objective, n.splits = 3, min.node.size = 1, improve.first.split = 0.1, improve.n.splits = 0.1){
  p = ncol(x)
  for(i in 1:p){
    effect = FeatureEffects$new(model, method = "ice", grid.size = 20, features = colnames(x)[i])
    # Get ICE values and arrange them in a horizontal matrix
    Y = spread(effect$results[[colnames(x)[i]]], .borders, .value)
    Y = Y[, setdiff(colnames(Y), c(".type", ".id", ".feature"))]
    for(j in 1:nrow(Y)){
      Y[j,] = as.numeric(unname(Y[j,])) - mean(as.numeric(unname(Y[j,])))
    }
    result.splits = sp(x = x, y = Y, feat = colnames(x)[i], n.splits = n.splits, min.node.size = min.node.size,
                       objective =objective, improve.first.split = improve.first.split, improve.n.splits = improve.n.splits)
    print(result.splits)
    
    if(!is.null(result.splits)){
      result.sub = data.frame(result.splits[[1]])
      result.sub$improvement = (result.splits[[2]] - result.sub$objective.value)/result.splits[[2]]
      result.sub$ICEfeature = colnames(x)[i]
      if(i ==1) result = result.sub
      else result = bind_rows(result, result.sub)
    }
    
  }
  return(result)
}



potential.interactions = find_potential_interactions(X, model, SS_fre, 3, min.node.size= 10, improve.first.split = 0.1, improve.n.splits = 0.1)
interactions = find_true_interactions(X, model)

find_true_interactions <- function(x, model, objective = SS_fre, n.splits = 3, min.node.size = 10, improve.first.split = 0.1, improve.n.splits = 0.1, interaction.factor = 0.05){
  potential.interactions = find_potential_interactions(x = x, model = model, objective = objective, n.splits = n.splits, min.node.size= min.node.size)
  
  
  true.interactions = potential.interactions[which(potential.interactions$improvement > interaction.factor),]
  
  ind = c()
  for(i in 1:nrow(true.interactions)){
    
    n.row = nrow(true.interactions[which(true.interactions$feature %in% true.interactions$ICEfeature[i]) %in% which(true.interactions$ICEfeature %in% true.interactions$feature[i]),])
    if(n.row == 0){
      ind = c(ind, i)
    }
  }
  true.interactions.all = true.interactions[-ind,]
  return(true.interactions.all[order(true.interactions.all$improvement, decreasing = TRUE),])
   
  
}


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




opt.feature = lapply(X[,-which(colnames(X)==feat)], function(xval) {
  optimizer(feat = feat, x = X, xval = xval, y = Y, n.splits = n.splits, min.node.size = min.node.size, 
            objective = objective)
})


