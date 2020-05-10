#' @title Performs a single split and measures the objective
#'
#' @description
#' Uses values in split.points as candidates to make intervals.
#' Sums up the objective in each interval (produced by the split.points).
#'
#' @param split.points [\code{vector}]\cr
#'   Either a single split point (for binary splits) or a vector of split points (for multiple splits)
#' @param xval feature values
#' @param y target
#' @return value of the objective
#' @export
#'test
perform_split = function(split.points, xval, y, min.node.size, objective) {
  UseMethod("perform_split")
}

# TODO: perform_split.numeric 
perform_split = function(split.points, xval, y, min.node.size, dist.between, lambda, objective) {
  
  
  # args = formalArgs(objective)
  # deparse(body(objective))
  # always increasing split points
  # split.points = xval[split.points] # integer optim
  split.points = sort.int(split.points)
  # assign intervalnr. according to split points
  node.number = findInterval(x = xval, split.points, rightmost.closed = TRUE) + 1
  # compute size of each childnode
  node.size = tabulate(node.number)
  # if minimum node size is violated, return Inf
  if (min(node.size) < min.node.size)
    return(Inf)
  # compute objective in each interval and sum it up
  y.list = split(y, node.number)
  # x.list only needed if this is used in the objective
  requires.x = formals(objective)[["requires.x"]]
  if (isTRUE(requires.x))
    x.list = split(xval, node.number) else
      x.list = NULL
  
  res = vapply(seq_along(y.list), FUN = function(i) {
    objective(y = y.list[[i]], x = x.list[[i]])
  }, FUN.VALUE = NA_real_, USE.NAMES = FALSE)
  
  # calculate cluster PDPs
  pdp_clust = sapply(seq_along(y.list), FUN = function(i) {
    colMeans(y.list[[i]])
  }, USE.NAMES = FALSE)
  
  
  # if regularization chosen to be "fre", then objective is adjusted accordingly (described in test_regularizationOfObjective.R)
  if (dist.between=="fre"){
    if (is.null(lambda)){
      lambda = 0.2
    }
    SS_betw = dist_fre(y=Y, x = xval, pdp = pdp_clust)
    return(sum(res) + lambda * sum(res)/(SS_betw/(length(split.points))))
  }
  # if regularization (dist.between) chosen to be "sd" adjustments of objective as described in test file
  else if (dist.between=="sd"){
    if (is.null(lambda)){
      lambda = 0.025
    }
    SD_dev = dist_sd(y=Y, x = x.list, pdp = pdp_clust)
    return(sum(res) + lambda * sum(res)/SD_dev)
  }
  
}


# distance measure for regularization method "fre"
dist_fre = function(y, x, pdp, requires.x = FALSE) { # slow
  # using only y-axis of curves is enough as x-axis is always the same for all curves
  require(kmlShape)
  center = colMeans(y)
  grid.x = as.numeric(names(center))
  pdp.y = unname(center)
  dist = apply(pdp, 2, function(ice) distFrechet(grid.x, pdp.y, grid.x, ice, FrechetSumOrMax = "sum"))
  sum(dist)
}

# distance measure for regularization method "sd"
dist_sd = function(y, x, pdp, requires.x = TRUE) { # slow
  # using only y-axis of curves is enough as x-axis is always the same for all curves
  require(kmlShape)
  center = colMeans(y)
  grid.x = as.numeric(names(center))
  pdp.y = unname(center)
  sd_sub = apply(pdp, 2, function(pdp_clust) sd(pdp_clust))
  weights = sapply(x, function(var) length(var))
  #weights = weights/sum(weights)
  abs(weighted.mean(sd_sub, w = unname(weights))-sd(pdp.y))
}

