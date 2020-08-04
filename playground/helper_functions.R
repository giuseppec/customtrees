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

# Frechet distance FDA measure
# SS_fre = function(y, x, requires.x = FALSE, ...) { # slow
#   # using only y-axis of curves is enough as x-axis is always the same for all curves
#   require(kmlShape)
#   center = colMeans(y)
#   grid.x = as.numeric(names(center))
#   pdp.y = unname(center)
#   dist = apply(y, 1, function(ice) distFrechet(grid.x, pdp.y, grid.x, ice, FrechetSumOrMax = "sum"))
#   sum(dist)
# }
# 
# # Frechet distance measure - with filtered ice curves
# SS_fre_filtered = function(y, x, sub.number, requires.x = FALSE, feat, x.all, ...) {
#   require(kmlShape)
#   # use only ice curves that are available for the combination of the two features -> no extrapolation
#   indices = filter(feat, x.all, y, sub.number)
#   y.filtered = y[,indices, drop = FALSE]
#   center = colMeans(y.filtered)
#   grid.x = as.numeric(names(center))
#   pdp.y = unname(center)
#   dist = apply(y.filtered, 1, function(ice) distFrechet(grid.x, pdp.y, grid.x, ice, FrechetSumOrMax = "sum"))
#   sum(dist)*20/length(indices)
# }

SS_fre = function(y, x, requires.x = FALSE, ...) { # slow
  require(Rfast)
  ypred = Rfast::colMedians(as.matrix(y))
  sum(t(abs(t(y) - ypred)))
}

# Frechet distance measure - with filtered ice curves
SS_fre_filtered = function(y, x, sub.number, requires.x = FALSE, feat, x.all, ...) {
  require(Rfast)
  ycols = ncol(y)
  # use only ice curves that are available for the combination of the two features -> no extrapolation
  indices = filter(feat, x.all, y, sub.number)
  y.filtered = y[,indices, drop = FALSE]
  dist = SS_fre(y.filtered)
  sum(dist)*ycols/length(indices)
}

# NEW:
# Filter function to use only ice curves within grid points to find best split point
# not yet included: categorical features (only numeric and one-hot-encoded)
# needs to be integrated in objective
filter = function(feat, x.all, Y, sub.number){
  values = unique(x.all[sub.number,feat])
  if (length(unique(x.all[,feat])) > 2) {
    grid.points = as.numeric(names(Y))
    break.points = grid.points[1:(length(grid.points) - 1)] + (grid.points[2:length(grid.points)] - grid.points[1:(length(grid.points) - 1)]) / 2
    range = cut(values, breaks = c(min(x.all[,feat]), break.points, max(x.all[,feat])), labels = c(names(Y)), include.lowest = TRUE, right = TRUE)
    return(which(names(Y) %in% unique(range)))
  }
  else if (length(unique(x.all[,feat])) == length(values)) {
    return(c(1:length(names(Y))))
  }
  else if (length(values) == 1 & length(unique(x.all[,feat])) == 2) {
    if (values < mean(unique(x.all[,feat]))) {
      return(c(1:round(ncol(Y)/2,0)))
    }
    else if (values > mean(unique(x.all[,feat]))) {
      return(c(round(ncol(Y)/2,0):ncol(Y)))
    }
  }
  
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
  if (extrapol == TRUE) {
    y.list.filtered = lapply(seq_along(y.list), FUN = function(i) {
      ind = filter(feat, X, Y, which(node.number == i))
      y.list[[i]][,ind, drop = FALSE]
    })
  }
  else y.list.filtered = y.list
  
  return(y.list.filtered)
}


# prepare data for plotting
plot.prep.ind = function(i, x){
  x$.id = 1:nrow(x)
  x = gather(x, .borders, .value, colnames(x)[1]:colnames(x)[ncol(x) - 1], factor_key = TRUE)
  x$.split = i
  return(x)
}

plot.prep.full = function(data){
  data.prep = lapply(seq_along(data), FUN = function(i) {
    plot.prep.ind(i, data[[i]])
  })
  data.prep = rbindlist(data.prep)
}

