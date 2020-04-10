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
#'
perform_split = function(split.points, xval, y, min.node.size, objective) {
  UseMethod("perform_split")
}

perform_split.numeric = function(split.points, xval, y, min.node.size, objective) {
  # always increasing split points
  # split.points = xval[split.points] # integer optim
  split.points = sort(split.points)
  # assign intervalnr. according to split points
  split.factor = findInterval(x = xval, split.points, rightmost.closed = TRUE)
  # compute size of each childnode
  node.size = table(split.factor)
  # if minimum node size is violated, return Inf
  if (any(node.size < min.node.size))
    return(Inf)
  # compute objective in each interval and sum it up
  y.list = split(y, split.factor)
  x.list = split(xval, split.factor)
  res = BBmisc::vnapply(seq_along(y.list), fun = function(i) {
    objective(x = x.list[[i]], y = y.list[[i]])
  })
  sum(res)
}

# not tested
perform_split.factor = function(split.points, xval, y, min.node.size, objective) {
  lev = levels(xval)
  xval = as.numeric(xval)
  split.points = which(lev %in% split.points)
  perform_split.numeric(xval = xval, split.points = split.points, y = y, min.node.size = min.node.size, objective = objective)
}

# perform_split2 = function(split.points, x, y, min.node.size = 10, objective) {
#   # always increasing split points
#   split.points = sort(split.points)
#   # assign intervalnr. according to split points
#   split.factor = findInterval(x, split.points, rightmost.closed = TRUE)
#   # compute size of each childnode
#   node.size = table(split.factor)
#   # if minimum node size is violated, return Inf
#   if (any(node.size < min.node.size))
#     return(Inf)
#   # compute objective in each interval and sum it up
#   d = data.table(x, y, split.factor)
#   sum(d[, .(obj = objective(x, y)), by = split.factor]$obj)
# }
