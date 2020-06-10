library(dplyr)
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
SS_fre = function(y, x, X, sub_number, feat, requires.x = FALSE) { # slow
  # using only y-axis of curves is enough as x-axis is always the same for all curves
  require(kmlShape)
  center = colMeans(y)
  grid.x = as.numeric(names(center))
  pdp.y = unname(center)
  dist = apply(y, 1, function(ice) distFrechet(grid.x, pdp.y, grid.x, ice, FrechetSumOrMax = "sum"))
  sum(dist)
}

# Frechet distance measure - with filtered ice curves
SS_fre_filtered = function(y, x, X, sub_number, requires.x = FALSE, feat) { 
  require(kmlShape)
  # use only ice curves that are available for the combination of the two features -> no extrapolation
  indices = filter(feat, X, y, sub_number)
  y.filtered = y[,indices, drop = FALSE]
  center = colMeans(y.filtered)
  grid.x = as.numeric(names(center))
  pdp.y = unname(center)
  dist = apply(y.filtered, 1, function(ice) distFrechet(grid.x, pdp.y, grid.x, ice, FrechetSumOrMax = "sum"))
  sum(dist)
}


# NEW:
# Filter function to use only ice curves within grid points to find best split point
# not yet included: categorical features (only numeric and one-hot-encoded)
# needs to be integrated in objective
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
      return(c(1:round(ncol(Y)/2,0)))
    }
    else if(values > mean(unique(X[,feat]))){
      return(c(round(ncol(Y)/2,0):ncol(Y)))
    }
  }
  
}



#------------------------------------------------------------------------------------------------------------
# functions for splitting


# Anpassung: "feat" (feature für welches ICE curves berechnet werden) muss beim Funktionsaufruf angegeben werden
# für filter Funktion in objective
# objective mit Filterung benötigt mehr Parameter - vereinfachung möglich?
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
  
  # x.list only needed if this is used in the objective
  requires.x = formals(objective)[["requires.x"]]
  if (isTRUE(requires.x))
    x.list = split(xval, node.number) else
      x.list = NULL
  
  res = vapply(seq_along(y.list), FUN = function(i) {
    objective(y = y.list[[i]], x = x.list[[i]], X = x, sub_number = which(node.number == i), feat = feat)
  }, FUN.VALUE = NA_real_, USE.NAMES = FALSE)
  sum(res)
}





# optimization function: test with binary split - return also ice curves for plotting
# Anpassung: 
# Parameter feat mit übergeben
# Parameter extrapol (True oder False, abhängig, ob Extrapolationsproblem beachtet werden soll oder nicht)
find_best_binary_split = function(feat, x, xval, y, n.splits = 1, min.node.size = 1, extrapol = TRUE,
                                  objective, ...) {
  assert_choice(n.splits, choices = 1)
  
  # use different split candidates to perform split
  q = generate_split_candidates(xval, use.quantiles = TRUE)
  splits = vapply(q, FUN = function(i) {
    perform_split(i, feat = feat, x = x, xval = xval, y = y, min.node.size = min.node.size, 
                  objective = objective)
  }, FUN.VALUE = NA_real_, USE.NAMES = FALSE)
  # select the split point yielding the minimal objective
  best = which.min(splits)
  # return also ice curves of best split point for plotting
  #ice = get_ice_curves(q[best], feat = feat, X = x, xval = xval, y = y, extrapol = extrapol)
  
  return(list(split.points = q[best], objective.value = splits[best]#, ice.curves = ice
              ))
}



# Optimization with MA-LS Chains: finds best split points and returns value of objective function
# Anpassung: 
# Parameter feat mit übergeben
# Parameter extrapol (True oder False, abhängig, ob Extrapolationsproblem beachtet werden soll oder nicht)
find_best_multiway_split2 = function(feat, x, xval, y, n.splits = 1, min.node.size = 1, extrapol = TRUE,
                                     objective, control = NULL, ...) {
  unique.x = length(unique(xval))
  # max. number of splits to be performed must be unique.x-1
  if (n.splits >= unique.x)
    n.splits = unique.x - 1
  # if only one split is needed, we use exhaustive search
  if (n.splits == 1)
    return(find_best_binary_split(feat, x, xval, y, n.splits = n.splits,
                                  min.node.size = min.node.size, extrapol = extrapol, objective = objective))
  
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
                  objective = objective)
  best = Rmalschains::malschains(.perform_split, lower = lower, upper = upper,
                                 initialpop = pop, verbosity = 0, control = control, ...) #   control = malschains.control(istep = 300, ls = "sw"),
  # return also ice curves of best split point for plotting
  #ice = get_ice_curves(sort(best$sol), feat = feat, X = x, xval = xval, y = y, extrapol = extrapol)
  
  return(list(split.points = sort(best$sol), objective.value = best$fitness#, ice.curves = ice
              ))
}


# Anpassung: 
# Parameter feat mit übergeben - für feat selbst wird optimizer nicht berechnet, da keine Interaktion zu sich selbst untersucht wird
# Parameter extrapol (True oder False, abhängig, ob Extrapolationsproblem beachtet werden soll oder nicht)
split_parent_node = function(Y, X, feat, n.splits = 1, min.node.size = 1, optimizer, extrapol = TRUE,
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
              extrapol = TRUE, objective = objective, ...) 
  })
  
  result = rbindlist(lapply(opt.feature, as.data.frame), idcol = "feature")
  result = result[, .(split.points = list(split.points)), by = c("feature", "objective.value"), with = TRUE]
  result$best.split = result$objective.value == min(result$objective.value)
  #result = result[, best.split := objective.value == min(objective.value)]
  return(result)
}

#------------------------------------------------------------------------------------------------------------
# functions for plotting


# get ice curves function for plotting
get_ice_curves <- function(Y, X, result, extrapol = TRUE){
  assert_data_table(result)
  # TODO: fix bug if more than one feature have the same best objective
  feature = unique(result$feature[result$best.split])
  split.points = unlist(result$split.points[result$best.split])
  split.points = sort.int(split.points)
  node.number = findInterval(x = X[,feature], split.points, rightmost.closed = TRUE) + 1
  
  # split y according to node.number
  y.list = split(Y, node.number)
  
  # ice curve feature
  feat = colnames(X)[which(!(colnames(X) %in% result$feature))]
  
  #filter ice curves in case extrapol = TRUE
  if(extrapol == TRUE){
    y.list.filtered = sapply(seq_along(y.list), FUN = function(i) {
      ind = filter(feat, X, Y, which(node.number == i))
      y.list[[i]][,ind, drop = FALSE]
    }, USE.NAMES = FALSE)
  }
  else y.list.filtered = y.list
  
  return(y.list.filtered)
}


# prepare data for plotting
plot.prep.ind = function(i, x){
  x$.id = 1:nrow(x)
  x = gather(x, .borders, .value, colnames(x)[1]:colnames(x)[ncol(x)-1], factor_key=TRUE)
  x$.split = i
  return(x)
}

plot.prep.full = function(data){
  data.prep = lapply(seq_along(data), FUN = function(i) {
    plot.prep.ind(i, data[[i]])
  })
  data.prep = rbindlist(data.prep)
}


#------------------------------------------------------------------------------------------------------------
# functions to find interactions 

# interaction effects for one feature
interactions_per_feature = function(x, y, feat, n.splits, min.node.size, objective.split, objective.total, improve.first.split, improve.n.splits){
  res = list()
  for(i in 1:n.splits) {
    res[[i]] = split_parent_node(Y = y, X = x, feat = feat, n.splits = i, min.node.size = min.node.size, 
                                 optimizer = find_best_multiway_split2, objective = objective.split)
    if(i == 1){
      SS_total = objective.total(y = y, x = x)
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


# find potential interactions between all features in X 
find_potential_interactions <- function(x, model,  objective.split, objective.total, n.splits = 3, min.node.size = 1, improve.first.split = 0.1, improve.n.splits = 0.1){
  p = ncol(x)
  for(i in 1:p){
    effect = FeatureEffects$new(model, method = "ice", grid.size = 20, features = colnames(x)[i])
    # Get ICE values and arrange them in a horizontal matrix
    Y = spread(effect$results[[colnames(x)[i]]], .borders, .value)
    Y = Y[, setdiff(colnames(Y), c(".type", ".id", ".feature"))]
    for(j in 1:nrow(Y)){
      Y[j,] = as.numeric(unname(Y[j,])) - mean(as.numeric(unname(Y[j,])))
    }
    result.splits = interactions_per_feature(x = x, y = Y, feat = colnames(x)[i], n.splits = n.splits, min.node.size = min.node.size,
                       objective.split =objective.split, objective.total = objective.total, 
                       improve.first.split = improve.first.split, improve.n.splits = improve.n.splits)
    
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


# find "true" interactions 
find_true_interactions <- function(x, model, objective.split = SS_fre_filtered, objective.total = SS_fre, n.splits = 3, min.node.size = 10, improve.first.split = 0.1, improve.n.splits = 0.1, interaction.factor = 0.05){
  # get potential interactions
  potential.interactions = find_potential_interactions(x = x, model = model, objective.split = objective.split, objective.total = objective.total, 
                                                       n.splits = n.splits, min.node.size= min.node.size)
  # user only interactions that have a higher effect than interaction.factor and that are two-sided
  true.interactions = potential.interactions[which(potential.interactions$improvement > interaction.factor),]
  
  ind = c()
  for(i in 1:nrow(true.interactions)){
    n.row = nrow(true.interactions[which(true.interactions$feature %in% true.interactions$ICEfeature[i]) %in% which(true.interactions$ICEfeature %in% true.interactions$feature[i]),])
    if(n.row == 0){
      ind = c(ind, i)
    }
  }
  if(length(ind) != 0) true.interactions = true.interactions[-ind,]
  
  return(true.interactions[order(true.interactions$improvement, decreasing = TRUE),])
}



#------------------------------------------------------------------------------------------------------------
# simulation example to test  functions


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
effect = FeatureEffects$new(model, method = "ice", grid.size = 20, features = "z")
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
# test plotting with extrapolation

res = split_parent_node(Y = Y, X = X, feat = "z", optimizer = find_best_multiway_split2, n.splits = 1, objective = SS_fre_filtered, min.node.size = 10)
ice = get_ice_curves(X = X, Y = Y, result = res)


plot.data = plot.prep.full(ice)

# Plot splitted ICE curves
ggplot(plot.data, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id)) + facet_grid(~ .split)



#------------------------------------------------------------------------------------------------------------
# test interaction functions


potential.interactions = find_potential_interactions(X, model, SS_fre_filtered, SS_fre, 3, min.node.size= 10, improve.first.split = 0.1, improve.n.splits = 0.1)
interactions = find_true_interactions(X, model)

# compare with hstatistics
ia1 = Interaction$new(model, feature = "z", grid.size = 20)
ia2 = Interaction$new(model, feature = "d", grid.size = 20)
ia3 = Interaction$new(model, feature = "x5", grid.size = 20)
ia4 = Interaction$new(model, feature = "x6", grid.size = 20)
ia = Interaction$new(model, grid.size = 20)




#-------------------------------------------------------------------------------------------------------------
# simulation studies with correlated features

source("R/helper_functions_sim.r")

# generate data
p = 10
n = 500

set.seed(1234)
Xsim = sim_data(N=n, n_groups= 1, y = 1, size_groups = 10, prop = c(0.4))
epsilon = rnorm(n)

Xsim = as.data.frame(Xsim)
colnames(Xsim) = paste0("V",1:p)

y = Xsim[,2] + 5*cos(Xsim[,3]*5)*Xsim[,6] + ifelse(Xsim[,10] <= mean(Xsim[,10]), I(8*Xsim[,3]),0) + epsilon


dat = data.frame(Xsim,y)
X = dat[, setdiff(colnames(dat), "y")]



# Fit model and compute ICE for z
mod = ranger(y ~ ., data = dat, num.trees = 500)
pred = predict.function = function(model, newdata) predict(model, newdata)$predictions
model = Predictor$new(mod, data = X, y = y, predict.function = pred)

# interactions # todo: Fehler bei "ind" nochmal testen
potential.interactions = find_potential_interactions(X, model, SS_fre_filtered, SS_fre, 3, min.node.size= 20, improve.first.split = 0.1, improve.n.splits = 0.1)
interactions = potential.interactions[order(potential.interactions$improvement, decreasing = TRUE),]
#saveRDS(interactions, "interactions.rds")
# zu hohe improvements durchweg - gefilterte objective funktioniert hier nicht - kann ungefilterte einfach 
# verwendet werden? nochmals durchdenken
# Interaktion zwischen x3 und x10 wird eindeutig gefunden. Nicht eindeutig für x3 und x6 (auch nicht bei HStatistik)
# Evtl. mit angepasster objective

# compare with hstatistics
ia = Interaction$new(model, feature = colnames(X)[i], grid.size = 20)




effect = FeatureEffects$new(model, method = "ice", grid.size = 20, features = "V6")
# Get ICE values and arrange them in a horizontal matrix
Y = spread(effect$results$V6, .borders, .value)
Y = Y[, setdiff(colnames(Y), c(".type", ".id", ".feature"))]
# centered ICE curves
for(i in 1:nrow(Y)){
  Y[i,] = as.numeric(unname(Y[i,])) - mean(as.numeric(unname(Y[i,])))
}


res = split_parent_node(Y, X, feat = "V6", n.splits = 1, min.node.size = 30, optimizer = find_best_multiway_split2, extrapol = TRUE, objective = SS_fre_filtered)

ice = get_ice_curves(X = X, Y = Y, result = res)


plot.data = plot.prep.full(ice)

# Plot splitted ICE curves
ggplot(plot.data, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id)) + facet_grid(~ .split)
