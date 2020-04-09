split_optimizer = function(x, y, n.splits = 1, min.node.size = 1, objective, ...) {
  init = rep(0, n.splits)

  # sample split points from observed x values
  gr = function(i, x, y, objective) {
    npar = length(i)
    q = quantile(x, seq(0.01, 0.99, by = 0.01), type = 1)
    par = sample(q, size = npar)
    sort(par)
  }

  # use simulated annealing to find optimal split points
  res = optim(init, fn = perform_split,
    x = x, y = y, objective = objective, min.node.size = min.node.size,
    method = "SANN", gr = gr, ...)

  list(split.points = res$par, objective.value = res$value)
}

split_optimizer_exhaustive = function(x, y, n.splits = 1, min.node.size = 1, objective, ...) {
  assert_choice(n.splits, choices = 1)

  q = quantile(x, seq(0.01, 0.99, by = 0.01), type = 1) # to speedup we use only quantile values as possible split points
  splits = vnapply(q, function(i) {
    perform_split(i, x = x, y = y, objective = objective, min.node.size = min.node.size)
  })
  # select the split point yielding the minimal objective
  best = which.min(splits)
  # split.index = which(x == q[best]) # TODO: potential bug if quantile value not observed value in x?
  list(split.points = q[best], objective.value = splits[best])
}

#
# opt = DEoptim(fn = perform_split, lower = rep(min(x), 1), upper = rep(max(x), 1),
#   x = x, y = y, objective = SS, control = DEoptim.control(initialpop = cbind(x)))
# opt$optim
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
