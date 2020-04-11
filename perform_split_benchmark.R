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
  yleft = y[left.ind] # TODO: error handling if vector is empty or contains too few values
  yright = y[!left.ind] # TODO: error handling if vector is empty or contains too few values

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
  sum(data[, .(obj = objective(xval, y)), by = node.number]$obj)
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
  d = data.table::data.table(x = xval, y = y, node.number = node.number)
  sum(d[, .(obj = objective(xval, y)), by = node.number]$obj)
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
