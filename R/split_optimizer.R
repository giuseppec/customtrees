# This is faster for binary splits (but does not work for multiple splits)
find_best_binary_split = function(xval, y, n.splits = 1, min.node.size = 10,
  objective, ...) {
  assert_choice(n.splits, choices = 1)

  # use different split candidates to perform split
  q = generate_split_candidates(xval, n.quantiles = 100, min.node.size = min.node.size)
  splits = vapply(q, FUN = function(i) {
    perform_split(i, xval = xval, y = y, min.node.size = min.node.size,
      objective = objective)
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
      min.node.size = min.node.size, objective = objective))

  # generate initial population
  xval.candidates = generate_split_candidates(xval, n.quantiles = NULL, min.node.size = min.node.size)
  # pop = replicate(n.splits, sample(xval.candidates, size = min(n.splits*10, length(xval.candidates))))
  pop.size = max(min(50, floor(length(xval.candidates)/n.splits)), 1)
  pop = t(replicate(pop.size, sort(sample(xval.candidates, size = n.splits))))
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
      objective = objective)
  best = Rmalschains::malschains(.perform_split, lower = lower, upper = upper,
    initialpop = pop, verbosity = 0, control = control, maxEvals = maxEvals,
    ...) #   control = malschains.control(istep = 300, ls = "sw"),

  return(list(split.points = sort(best$sol), objective.value = best$fitness))
}

# Optimization with Hooke-Jeeves and restarts
find_best_multiway_split_hjn = function(xval, y, n.splits = 1, min.node.size = 10,
  objective, control = NULL, ...) {
  n.splits = adjust_nsplits(xval, n.splits)
  # if only one split is needed, we use exhaustive search
  if (n.splits == 1)
    return(find_best_binary_split(xval, y, n.splits = n.splits, min.node.size = min.node.size,
      objective = objective))

  if (is.null(control))
    control = list(maxfeval = 640, maxrestarts = 3)

  # find optimal split points constrained optimization
  n.quantiles = NULL #min(4*n.splits, 40)
  candidate = generate_split_candidates(xval, n.quantiles = n.quantiles, min.node.size = min.node.size)

  lower = rep(min(xval), n.splits)
  upper = rep(max(xval), n.splits)
  #init = lhs::randomLHS(control$maxrestarts, n.splits)
  #init = apply(init, 2, sort)
  #init = apply(init, 1, function(x) lower[1] + (upper[1] - lower[1])*x)
  init = replicate(control$maxrestarts, sort(sample(candidate, size = n.splits)))
  #init = sort(sample(candidate, size = n.splits))

  best = vector("list", length = control$maxrestarts)
  i = 1
  funeval = 0
  while (i <= control$maxrestarts & funeval < control$maxfeval) {
    best[[i]] = optimx::hjn(par = init[, i], fn = perform_split, lower = lower, upper = upper,
      xval = xval, y = y, min.node.size = min.node.size, objective = objective, control = control)
    funeval = funeval + best[[i]]$counts[1]
    i = i + 1
  }
  best = best[[which.min(unlist(lapply(best, function(x) x$value)))]]

  return(list(split.points = sort(best$par), objective.value = best$value))
}

# Optimization with Hooke-Jeeves derivative-free (fast, for low dimensions gets stuck in local minima -> restarts?)
find_best_multiway_split_hjkb = function(xval, y, n.splits = 1, min.node.size = 10,
  objective, control = NULL, ...) {
  n.splits = adjust_nsplits(xval, n.splits)
  # if only one split is needed, we use exhaustive search
  if (n.splits == 1)
    return(find_best_binary_split(xval, y, n.splits = n.splits, min.node.size = min.node.size,
      objective = objective))

  if (is.null(control))
    control = list(maxfeval = 640, tol = 1e-4)

  # find optimal split points constrained optimization
  n.quantiles = NULL #min(4*n.splits, 40)
  candidate = generate_split_candidates(xval, n.quantiles = n.quantiles, min.node.size = min.node.size)
  init = sort(sample(candidate, size = n.splits))
  # init = lapply(1:10, function(i) sort(sample(candidate, size = n.splits))) #candidate[1:n.splits]
  # init.eval = vapply(init, perform_split,
  #   xval = xval, y = y, min.node.size = min.node.size, objective = objective,
  #   FUN.VALUE = NA_real_, USE.NAMES = FALSE)
  # init = init[[which.min(init.eval)]]

  lower = rep(min(xval), n.splits)
  upper = rep(max(xval), n.splits)
  best = dfoptim::hjkb(par = init, fn = perform_split, lower = lower, upper = upper, control = control,
    xval = xval, y = y, min.node.size = min.node.size, objective = objective)

  # best = dfoptim::hjkb(par = sort(best$par), fn = perform_split, lower = lower, upper = upper, control = control,
  #   xval = xval, y = y, min.node.size = min.node.size, objective = objective)

  return(list(split.points = sort(best$par), objective.value = best$value))
}

# Optimization with GenSA
find_best_multiway_split_gensa = function(xval, y, n.splits = 1, min.node.size = 10,
  objective, control = NULL, ...) {
  n.splits = adjust_nsplits(xval, n.splits)
  # if only one split is needed, we use exhaustive search
  if (n.splits == 1)
    return(find_best_binary_split(xval, y, n.splits = n.splits,
      min.node.size = min.node.size, objective = objective))

  # generate initial population
  n.quantiles = NULL
  candidate = generate_split_candidates(xval, n.quantiles = n.quantiles, min.node.size = min.node.size)
  init = sort(sample(candidate, size = n.splits))

  if (is.null(control))
    control = list(max.call = 640, nb.stop.improvement = min(n.splits*10, 100), smooth = FALSE)

  # possible lower and upper values
  lower = rep(min(xval), n.splits)
  upper = rep(max(xval), n.splits)

  # optimization function
  best = GenSA::GenSA(par = init, fn = perform_split, lower = lower, upper = upper,
    control = control, xval = xval, y = y, min.node.size = min.node.size, objective = objective)

  return(list(split.points = sort(best$par), objective.value = best$value))
}

# helper functions
adjust_nsplits = function(xval, n.splits) {
  # max. number of splits to be performed must be unique.x-1
  unique.x = length(unique(xval))
  if (n.splits >= unique.x)
    n.splits = unique.x - 1
  return(n.splits)
}

get_closest_point = function(split.points, xval) {
  # get closest value from xval
  xval = sort(xval)
  # xval = xval[-c(1, length(xval))]
  ind = vapply(split.points, function(i) {
    d = i - xval
    which.min(abs(d))
    #which.min(d[d >= 0])
  }, FUN.VALUE = NA_real_, USE.NAMES = FALSE)
  return(xval[ind])
}

adjust_split_point = function(split.points, xval) {
  # use a value between two subsequent points
  q = split.points
  x.unique = sort(unique(xval))
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
  xval = sort(xval)
  # try to ensure min.node.size between points (is not guaranteed)
  chunk.ind = seq.int(min.node.size + 1, length(xval) - min.node.size, by = min.node.size)
  xadj = unique(quantile(xval, prob = chunk.ind/length(xval), type = 1))
  # xadj = xval[chunk.ind]

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

# # finds best split points and returns value of objective function
# find_best_multiway_split2 = function(xval, y, n.splits = 1, min.node.size = 10 objective, use.quantiles = TRUE, ...) {
#   # use simulated annealing to find optimal split points
#   init = rep(0, n.splits)
#   nlm(f = perform_split, p = init, xval = xval, y = y, objective = objective, min.node.size = min.node.size)
#   fun = function(s) perform_split(split.points = s, xval = x, y = y, objective = objective, min.node.size = min.node.size)
#   best = optimization::optim_sa(fun = fun, start = init, lower = rep(min(x), n.splits), upper = rep(max(x), n.splits))
#   return(list(split.points = best$par, objective.value = best$function_value))
# }
#
# opt = DEoptim(fn = perform_split, lower = rep(min(x), 1), upper = rep(max(x), 1), min.node.size = 100,
#   x = x, y = y, objective = SS, control = DEoptim.control(initialpop = cbind(x)))
# opt$optim
#
# MU = 99L; LAMBDA = 5L; MAX.ITER = 200L
# opt = ecr(fitness.fun = perform_split, min.node.size = 100, x = x, y = y, objective = SS,
#   minimize = TRUE,
#   mutator = setup(mutInsertion, ind = 1:length(x)),
#   n.objectives = 1, mu = MU, lambda = LAMBDA,
#   representation = "float"#,
#   #initial.solutions = as.list(quantile(x, 1:99/100, type = 1)),
#   #custom.constants = as.list(x)
#   )
# opt
# opt$best.x
#
# MU = 25; LAMBDA = 30; MAX.ITER = 100
#
# ctrl = initECRControl(fitness.fun = perform_split, n.objectives = 1L, minimize = TRUE)
# ctrl = registerECROperator(ctrl, "selectForSurvival", selGreedy)
# population = gen(expr = {sample(x, size = 1)}, n = 100)
# fitness = evaluateFitness(ctrl, population, min.node.size = 10, x = x, y = y, objective = SS)
#
# for (i in seq_len(MAX.ITER)) {
#   offspring = generateOffspring(ctrl, population, fitness, lambda = LAMBDA, p.recomb = 0.7, p.mut = 0.3)
#   fitness.o = evaluateFitness(ctrl, offspring)
#   sel = replaceMuPlusLambda(ctrl, population, offspring, fitness, fitness.o)
#   population = sel$population
#   fitness = sel$fitness
# }
# print(population[[which.max(fitness)]])
# print(max(fitness))
#
# pointReplace = makeMutator(
#   mutator = function(ind) {
#     idx = which(runif(nrow(ind)) < 0.1)
#     ind[idx ,] = matrix(runif(2*length(idx)), ncol = 2)
#     return(ind)
#     }, supported = custom)
#
#
# opt2 = find_best_multiway_split(x, y, objective = SS, n.splits = 2)
# opt2
#
# k = 1
# g = GenSA(rep(0, k), fn = perform_split, lower = rep(min(x), k), upper = rep(max(x), k), x = x, y = y, objective = SS)
#
# ga = ga(type = "real-valued", fitness = perform_split,
#   lower = rep(min(x), k), upper = rep(max(x), k), x = x, y = y, objective = SS, optim = TRUE,
#   optimArgs = list(method = "SANN"))
# ga@population
#
# library(NMOF)
# ?LSopt
#
# xval = x
#
# neighbour = function(x, ...) {
#   par = list(...)
#   sample(x, size = length(par[[1]]))
# }
#
# neighbour <- function(x, ...)
#   x + sample(seq(-3,3), length(x), replace=TRUE)
#
# LSopt(OF = perform_split, algo = list(x0 = rep(0, 1), neighbour = get_neighbour_function),
#   x = x, y = y, objective = SS)$xbest
#
