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
# perform_split = function(split.points, xval, y, min.node.size, objective) {
#   UseMethod("perform_split")
# }


# Anpassung: "feat" (feature für welches ICE curves berechnet werden) muss beim Funktionsaufruf angegeben werden
# für filter Funktion in objective
# objective mit Filterung benötigt mehr Parameter - vereinfachung möglich?
perform_split = function(split.points, feat, x, xval, y, min.node.size, objective,...) {
  
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
    objective(y = y.list[[i]], xval = x.list[[i]], x = x, sub.number = which(node.number == i), feat = feat,...)
  }, FUN.VALUE = NA_real_, USE.NAMES = FALSE)
  sum(res)
}



