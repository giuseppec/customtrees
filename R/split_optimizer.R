# finds best split points and returns value of objective function
find_best_multiway_split = function(xval, y, n.splits = 1, min.node.size = 1, objective, use.quantiles = TRUE, ...) {
  # sample split points from observed x values as candidates in optimization procedure
  gr = function(i, xval, y, objective, min.node.size) {
    npar = length(i)
    q = generate_split_candidates(xval, use.quantiles = use.quantiles)
    par = sample(q, size = npar)
    sort(par)
  }

  # use simulated annealing to find optimal split points
  init = rep(0, n.splits)
  best = optim(init, fn = perform_split,
    xval = xval, y = y, objective = objective, min.node.size = min.node.size,
    method = "SANN", gr = gr, ...)

  return(list(split.points = best$par, objective.value = best$value))
}

# This is faster for binary splits (but does not work for multiple splits)
find_best_binary_split = function(xval, y, n.splits = 1, min.node.size = 1, objective, use.quantiles = TRUE, ...) {
  assert_choice(n.splits, choices = 1)

  # use different split candidates to perform split
  q = generate_split_candidates(xval, use.quantiles = use.quantiles)
  splits = BBmisc::vnapply(q, function(i) {
    perform_split(i, xval = xval, y = y, objective = objective, min.node.size = min.node.size)
  })
  # select the split point yielding the minimal objective
  best = which.min(splits)

  return(list(split.points = q[best], objective.value = splits[best]))
}

generate_split_candidates = function(xval, use.quantiles = TRUE) {
  if (use.quantiles) { # to speedup we use only quantile values as possible split points
    q = unique(quantile(xval, seq(0.01, 0.99, by = 0.01), type = 1))
  } else {
    q = sort(unique(xval))
  }
  # use mid-point of two subsequent split candidates
  q = q + c(diff(q)/2, 0)
  return(q)
}

#
# opt = DEoptim(fn = perform_split, lower = rep(min(x), 1), upper = rep(max(x), 1), min.node.size = 10,
#   x = x, y = y, objective = SS, control = DEoptim.control(initialpop = cbind(x)))
# opt$optim
#
# MU = 99L; LAMBDA = 5L; MAX.ITER = 200L
# opt = ecr(fitness.fun = perform_split, min.node.size = 10, x = x, y = y, objective = SS,
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
