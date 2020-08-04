# This is faster for binary splits (but does not work for multiple splits)
find_best_binary_split = function(xval, y, n.splits = 1, min.node.size = 10,
  objective, ...) {
  assert_choice(n.splits, choices = 1)

  # use different split candidates to perform split
  q = generate_split_candidates(xval, n.quantiles = 100, min.node.size = min.node.size)
  splits = vapply(q, FUN = function(i) {
    perform_split(i, xval = xval, y = y, min.node.size = min.node.size,
      objective = objective, ...)
  }, FUN.VALUE = NA_real_, USE.NAMES = FALSE)
  # select the split point yielding the minimal objective
  best = which.min(splits)

  return(list(split.points = q[best], objective.value = splits[best]))
}

# finds best split points and returns value of objective function
# choose fastest according to https://www.researchgate.net/publication/311445822_Memetic_Algorithms_with_Local_Search_Chains_in_R_The_Rmalschains_Package
find_best_multiway_split = function(xval, y, n.splits = 1, min.node.size = 10,
  objective, control = NULL, ...) {
  find_best_multiway_split_mals(xval, y, n.splits, min.node.size, objective, control, ...)
}

# Optimization with MA-LS Chains
find_best_multiway_split_mals = function(xval, y, n.splits = 1, min.node.size = 10,
  objective, control = NULL, ...) {
  n.splits = adjust_nsplits(xval, n.splits)
  # if only one split is needed, we use exhaustive search
  if (n.splits == 1)
    return(find_best_binary_split(xval, y, n.splits = n.splits,
      min.node.size = min.node.size, objective = objective, ...))

  # generate initial population
  xval.candidates = generate_split_candidates(xval, n.quantiles = NULL, min.node.size = min.node.size)
  # pop = replicate(n.splits, sample(xval.candidates, size = min(n.splits*10, length(xval.candidates))))
  pop.size = max(min(50, floor(length(xval.candidates)/n.splits)), 1)
  pop = t(replicate(pop.size, sort.int(sample(xval.candidates, size = n.splits))))
  pop = unique(pop)

  if (is.null(control))
    control = Rmalschains::malschains.control(istep = 80, ls = "sw", threshold = 1e-4)
  maxEvals = 8*control$istep
  # # replace maximum number of evaluations depending of number of unique x values
  # xval.ncomb = prod(rep(length(xval.candidates), n.splits))
  # control$istep = min(control$istep, xval.ncomb)

  # possible lower and upper values
  lower = rep(min(xval.candidates), n.splits)
  upper = rep(max(xval.candidates), n.splits)

  # optimization function
  .perform_split = function(par)
    perform_split(split.points = par, xval = xval, y = y, min.node.size = min.node.size,
      objective = objective, ...)
  best = Rmalschains::malschains(.perform_split, lower = lower, upper = upper,
    initialpop = pop, verbosity = 0, control = control, maxEvals = maxEvals) #   control = malschains.control(istep = 300, ls = "sw"),

  return(list(split.points = sort.int(best$sol), objective.value = best$fitness))
}

# helper functions
adjust_nsplits = function(xval, n.splits) {
  # max. number of splits to be performed must be unique.x-1
  unique.x = length(unique(xval))
  if (n.splits >= unique.x)
    n.splits = unique.x - 1
  return(n.splits)
}

# get_closest_point2 = function(split.points, xval, min.node.size = 10) {
#   p = floor(length(xval)/min.node.size)
#   ret = unique(quantile(xval, probs = (1:(p - 1))/p))
#   vapply(ret, function(i) xval[which.min(abs(xval - i))], FUN.VALUE = NA_real_, USE.NAMES = FALSE)
# }

# replace split.points with closest value from xval taking into account min.node.size
get_closest_point = function(split.points, xval, min.node.size = 10) {
  xval = sort.int(xval)
  # try to ensure min.node.size between points (is not guaranteed if many duplicated values exist)
  chunk.ind = seq.int(min.node.size + 1, length(xval) - min.node.size, by = min.node.size)
  xadj = unique(xval[chunk.ind]) # unique(quantile(xval, prob = chunk.ind/length(xval), type = 1))
  # xval = xval[-c(1, length(xval))]
  split.adj = numeric(length(split.points))
  for (i in seq_along(split.adj)) {
    d = xadj - split.points[i]
    ind.closest = which.min(abs(d))
    split.adj[i] = xadj[ind.closest]
    xadj = xadj[-ind.closest] # remove already chosen value
  }
  # ind = vapply(split.points, function(i) {
  #   d = xadj - i
  #   res = which.min(abs(d))
  #   #which.min(d[d >= 0])
  #   return(res)
  # }, FUN.VALUE = NA_real_, USE.NAMES = FALSE)
  return(sort.int(split.adj))
}

adjust_split_point = function(split.points, xval) {
  # use a value between two subsequent points
  q = split.points
  x.unique = sort.int(unique(xval))
  ind = which(x.unique %in% q)
  ind = ind[ind < length(x.unique)]
  if (length(ind) != length(q)) {
    eps = min(diff(x.unique))/2
  } else {
    eps = (x.unique[ind + 1] - x.unique[ind])/2
  }
  q = q + eps #+ c(diff(q)/2, 0)
  q[q < x.unique[2]] = mean(x.unique[1:2])
  q[q > x.unique[length(x.unique) - 1]] = mean(x.unique[(length(x.unique) - 1):length(x.unique)])
  #q = q[q <= max(xval) & q >= min(xval)]
  return(q)
}

generate_split_candidates = function(xval, n.quantiles = NULL, min.node.size = 10) {
  assert_integerish(min.node.size, upper = floor((length(xval) - 1)/2))
  xval = sort.int(xval)
  # try to ensure min.node.size between points (is not guaranteed)
  chunk.ind = seq.int(min.node.size + 1, length(xval) - min.node.size, by = min.node.size)
  #xadj = unique(quantile(xval, prob = chunk.ind/length(xval), type = 1))
  xadj = xval[chunk.ind]

  if (!is.null(n.quantiles)) {
    # to speedup we use only quantile values as possible split points
    # qprobs = seq(1/n.quantiles, (n.quantiles - 1)/n.quantiles, by = 1/n.quantiles)
    qprobs = seq(0, 1, by = 1/n.quantiles)
    q = unique(quantile(xadj, qprobs, type = 1))
  } else {
    q = unique(xadj)
  }

  # use a value between two subsequent points
  q = adjust_split_point(q, xval)
  # x.unique = unique(xval)
  # ind = which(x.unique %in% q)
  # ind = ind[ind < length(x.unique)]
  # if (length(ind) != length(q)) {
  #   eps = min(diff(x.unique))/2
  # } else {
  #   eps = (x.unique[ind + 1] - x.unique[ind])/2
  # }
  # q = q + eps #+ c(diff(q)/2, 0)
  # q = q[q <= max(xval) & q >= min(xval)]
  return(q)
}
