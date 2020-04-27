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
#   node.number = findInterval(x, split.points, rightmost.closed = TRUE)
#   # compute size of each childnode
#   node.size = table(node.number)
#   # if minimum node size is violated, return Inf
#   if (any(node.size < min.node.size))
#     return(Inf)
#   # compute objective in each interval and sum it up
#   d = data.table(x, y, node.number)
#   sum(d[, .(obj = objective(x, y)), by = node.number]$obj)
# }
