
# cumsum representation
find_best_multiway_split_mals = function(xval, y, n.splits = 1, min.node.size = 10,
  objective, control = NULL, ...) {
  n.splits = adjust_nsplits(xval, n.splits)
  # if only one split is needed, we use exhaustive search
  if (n.splits == 1)
    return(find_best_binary_split(xval, y, n.splits = n.splits,
      min.node.size = min.node.size, objective = objective))

  if (is.null(control))
    control = Rmalschains::malschains.control(istep = 100, ls = "sw", threshold = 1e-4)
  maxEvals = 10*control$istep
  # # replace maximum number of evaluations depending of number of unique x values
  # xval.ncomb = prod(rep(length(xval.candidates), n.splits))
  # control$istep = min(control$istep, xval.ncomb)

  # possible lower and upper values
  lower = c(min(xval), rep(min(diff(sort(unique(xval)))), n.splits - 1))
  upper = c(max(xval), rep(diff(range(xval)), n.splits - 1))

  # optimization function
  .perform_split = function(par)
    perform_split(split.points = par, xval = xval, y = y, min.node.size = min.node.size,
      objective = objective)
  best = Rmalschains::malschains(.perform_split, lower = lower, upper = upper,
    verbosity = 0, control = control, maxEvals = maxEvals,
    ...) #   control = malschains.control(istep = 300, ls = "sw"),

  return(list(split.points = cumsum(best$sol), objective.value = best$fitness))
} # cumsum representation
perform_split = function(split.points, xval, y, min.node.size, objective) {
  # args = formalArgs(objective)
  # deparse(body(objective))
  # always increasing split points
  # split.points = xval[split.points] # integer optim

  split.points = cumsum(split.points)
  # assign intervalnr. according to split points

  node.number = findInterval(x = xval, split.points, rightmost.closed = TRUE) + 1
  # compute size of each childnode
  node.size = tabulate(node.number)
  # if minimum node size is violated, return Inf
  # TODO: instead of returning Inf try to avoid that this happens by fixing split points
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





SS_fre2 = function(y, x, requires.x = FALSE) { # slow
  # using only y-axis of curves is enough as x-axis is always the same for all curves
  require(kmlShape)
  center = colMeans(y)
  grid.x = as.numeric(names(center))
  pdp.y = unname(center)
  sum(y[, .(dist = distFrechet(grid.x, pdp.y, grid.x, .SD)), by = 1:nrow(y)]$dist)
}

perform_split0 = function(split.points, xval, y, min.node.size, objective) {
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
  x.list = split(xval, node.number)
  res = BBmisc::vnapply(seq_along(y.list), fun = function(i) {
    objective(x = x.list[[i]], y = y.list[[i]])
  })
  sum(res)
}

# works only for binary splits
perform_split1 = function(split.points, xval, y, min.node.size = 10, objective) {
  left.ind = xval < split.points
  # split y space
  if (is.data.frame(y)) {
    yleft = y[left.ind, ] # TODO: error handling if vector is empty or contains too few values
    yright = y[!left.ind, ] # TODO: error handling if vector is empty or contains too few values
  } else {
    yleft = y[left.ind] # TODO: error handling if vector is empty or contains too few values
    yright = y[!left.ind] # TODO: error handling if vector is empty or contains too few values
  }

  # split xval space (not required if objective is independend from xval)
  xleft = xval[left.ind]
  xright = xval[!left.ind]

  objective(x = xleft, y = yleft) + objective(x = xright, y = yright)
}

# only works univariate targets and features?
perform_split2 = function(split.points, data, feature, target, min.node.size = 10, objective) {
  y = data[[target]]
  xval = data[[feature]]
  # always increasing split points
  split.points = sort(split.points)
  # assign intervalnr. according to split points
  node.number = findInterval(x = xval, split.points, rightmost.closed = TRUE) + 1
  # compute size of each childnode
  node.size = tabulate(node.number)
  # if minimum node size is violated, return Inf
  if (min(node.size) < min.node.size)
    return(Inf)
  # compute objective in each interval and sum it up
  sum(data[, .(obj = objective(y = y, x = x)), by = node.number]$obj)
}

# This is faster for binary splits (but does not work for multiple splits)
find_best_binary_split2 = function(xval, y, n.splits = 1, min.node.size = 1,
  objective, ...) {
  assert_choice(n.splits, choices = 1)

  # use different split candidates to perform split
  q = generate_split_candidates(xval, n.quantiles = 100, min.node.size = min.node.size)
  splits = BBmisc::vnapply(q, function(i) {
    perform_split3(i, xval = xval, y = y, min.node.size = min.node.size,
      objective = objective)
  })
  # select the split point yielding the minimal objective
  best = which.min(splits)

  return(list(split.points = q[best], objective.value = splits[best]))
}

perform_split3 = function(split.points, xval, y, min.node.size = 10, objective) {
  # always increasing split points
  split.points = sort(split.points)
  # assign intervalnr. according to split points
  node.number = findInterval(x = xval, split.points, rightmost.closed = TRUE) + 1
  # compute size of each childnode
  node.size = tabulate(node.number)
  # if minimum node size is violated, return Inf
  if (min(node.size) < min.node.size)
    return(Inf)
  # compute objective in each interval and sum it up
  d = data.table::data.table(x = xval, y = list(y), node.number = node.number)
  sum(d[, .(obj = objective(y = y, x = x)), by = node.number]$obj)
}
nsim = 10000
set.seed(31415)
x = sort(runif(n = nsim, min = 0, max = 5*pi))
y = 4 + 2 * cos(x) + rnorm(nsim, mean = 0, sd = abs(cos(x)) / 2)
SS = function(x, y) {
  ypred = mean(y)
  sum((y - ypred)^2)
}

q = quantile(x, seq(0, 1, length.out = 100), type = 1)

dat = data.table::data.table(x = x, y = y)


library(microbenchmark)
microbenchmark(
  perform_split.numeric(q[3], x, y, objective = SS, min.node.size = 1),
  perform_split0(q[3], x, y, objective = SS, min.node.size = 1),
  perform_split1(q[3], x, y, objective = SS, min.node.size = 1),
  perform_split2(q[3], data = dat, target = "y", feature = "x", objective = SS, min.node.size = 1),
  perform_split3(q[3], x, y, objective = SS, min.node.size = 1), times = 100
)


cuts = sort(sample(q, size = 2))
microbenchmark(
  perform_split.numeric(cuts, x, y, objective = SS, min.node.size = 1),
  perform_split0(cuts, x, y, objective = SS, min.node.size = 1),
  perform_split2(cuts, data = dat, target = "y", feature = "x", objective = SS, min.node.size = 1),
  perform_split3(cuts, x, y, objective = SS, min.node.size = 1), times = 100
)

cat = findInterval(x, cuts, rightmost.closed = TRUE, all.inside = TRUE)
microbenchmark(
  Table(cat, names = FALSE), tabulate(cat), times = 100
)

node.size = tabulate(cat)
min.node.size = 3
microbenchmark(
  any(node.size < min.node.size),
  min(node.size) < min.node.size,
  node.size[which.min(node.size)] < min.node.size
)


library(profvis)
profvis(
  perform_split.numeric(q[c(3,10,18)], x, y, objective = SS, min.node.size = 1)
  )
