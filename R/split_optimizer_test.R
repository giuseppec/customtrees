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
  init = replicate(control$maxrestarts, sort.int(sample(candidate, size = n.splits)))
  #init = sort.int(sample(candidate, size = n.splits))
  
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
  
  return(list(split.points = sort.int(best$par), objective.value = best$value))
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
  init = sort.int(sample(candidate, size = n.splits))
  # init = lapply(1:10, function(i) sort.int(sample(candidate, size = n.splits))) #candidate[1:n.splits]
  # init.eval = vapply(init, perform_split,
  #   xval = xval, y = y, min.node.size = min.node.size, objective = objective,
  #   FUN.VALUE = NA_real_, USE.NAMES = FALSE)
  # init = init[[which.min(init.eval)]]
  
  lower = rep(min(xval), n.splits)
  upper = rep(max(xval), n.splits)
  best = dfoptim::hjkb(par = init, fn = perform_split, lower = lower, upper = upper, control = control,
    xval = xval, y = y, min.node.size = min.node.size, objective = objective)
  
  # best = dfoptim::hjkb(par = sort.int(best$par), fn = perform_split, lower = lower, upper = upper, control = control,
  #   xval = xval, y = y, min.node.size = min.node.size, objective = objective)
  
  return(list(split.points = sort.int(best$par), objective.value = best$value))
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
  init = sort.int(sample(candidate, size = n.splits))
  
  if (is.null(control))
    control = list(max.call = 640, nb.stop.improvement = min(20*n.splits, 100), smooth = FALSE)
  
  # possible lower and upper values
  lower = rep(min(xval), n.splits)
  upper = rep(max(xval), n.splits)
  
  # optimization function
  best = GenSA::GenSA(par = init, fn = perform_split, lower = lower, upper = upper,
    control = control, xval = xval, y = y, min.node.size = min.node.size, objective = objective)
  
  return(list(split.points = sort.int(best$par), objective.value = best$value))
}


# Optimization with
find_best_multiway_split_mbo = function(xval, y, n.splits = 1, min.node.size = 10,
  objective, control = NULL, ...) {
  n.splits = adjust_nsplits(xval, n.splits)
  # if only one split is needed, we use exhaustive search
  if (n.splits == 1)
    return(find_best_binary_split(xval, y, n.splits = n.splits, min.node.size = min.node.size,
      objective = objective))

  # find optimal split points constrained optimization
  n.quantiles = NULL #min(4*n.splits, 40)
  candidate = generate_split_candidates(xval, n.quantiles = n.quantiles, min.node.size = min.node.size)
  init = sort.int(sample(candidate, size = n.splits))


  obj.fun = makeSingleObjectiveFunction(
    name = "perform_split",
    fn = function(x) perform_split(sort.int(x), xval = xval, y = y, min.node.size = min.node.size, objective = objective),
    par.set = makeNumericParamSet(lower = min(xval), upper = max(xval), len = n.splits),
    global.opt.value = 0
  )

  ctrl = makeMBOControl(propose.points = 1)
  ctrl = setMBOControlTermination(ctrl, iters = 50L)
  ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI(),
    opt = "focussearch", opt.focussearch.points = 100L)

  lrn = makeMBOLearner(ctrl, obj.fun)
  lrn = setHyperPars(lrn, control = list(trace = FALSE))

  design = generateDesign(4*n.splits, getParamSet(obj.fun), fun = lhs::maximinLHS)

  best = mbo(obj.fun, design = design, learner = lrn, control = ctrl, show.info = T)

  return(list(split.points = sort.int(best$x), objective.value = best$y))
}



# Optimization with
find_best_multiway_split_hj = function(xval, y, n.splits = 1, min.node.size = 10,
  objective, control = NULL, ...) {
  n.splits = adjust_nsplits(xval, n.splits)
  # if only one split is needed, we use exhaustive search
  if (n.splits == 1)
    return(find_best_binary_split(xval, y, n.splits = n.splits, min.node.size = min.node.size,
      objective = objective))

  # find optimal split points constrained optimization
  n.quantiles = NULL #min(4*n.splits, 40)
  candidate = generate_split_candidates(xval, n.quantiles = n.quantiles, min.node.size = min.node.size)
  init = sort.int(sample(candidate, size = n.splits))

  lower = rep(min(xval), n.splits)
  upper = rep(max(xval), n.splits)

  best = adagio::hookejeeves(x0 = init, f = perform_split,
    xval = xval, y = y, min.node.size = min.node.size, objective = objective,
    lb = lower, ub = upper, maxfeval = 640, tol = 1e-4)

  return(list(split.points = sort.int(best$xmin), objective.value = best$fmin))
}

# Optimization with
find_best_multiway_split_ea = function(xval, y, n.splits = 1, min.node.size = 10,
  objective, control = NULL, ...) {
  n.splits = adjust_nsplits(xval, n.splits)
  # if only one split is needed, we use exhaustive search
  if (n.splits == 1)
    return(find_best_binary_split(xval, y, n.splits = n.splits, min.node.size = min.node.size,
      objective = objective))

  # find optimal split points constrained optimization
  n.quantiles = NULL #min(4*n.splits, 40)
  candidate = generate_split_candidates(xval, n.quantiles = n.quantiles, min.node.size = min.node.size)
  init = sort.int(sample(candidate, size = n.splits))

  lower = rep(min(xval), n.splits)
  upper = rep(max(xval), n.splits)

  # best = adagio::simpleEA(fn = perform_split,
  #   xval = xval, y = y, min.node.size = min.node.size, objective = objective,
  #   lower = lower, upper = upper, N = 50, confined = TRUE)
  best = adagio::pureCMAES(par = init, fun = perform_split,
    xval = xval, y = y, min.node.size = min.node.size, objective = objective,
    lower = lower, upper = upper)

  return(list(split.points = sort.int(best$xmin), objective.value = best$fmin))
}

# Optimization with
find_best_multiway_split_deopt = function(xval, y, n.splits = 1, min.node.size = 10,
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
  init = sort.int(sample(candidate, size = n.splits))

  lower = rep(min(xval), n.splits)
  upper = rep(max(xval), n.splits)

  repair = function(par, ...) {
    li = list(...)
    res = adjust_split_point(get_closest_point(par, li$xval), li$xval)
    #cat(res, fill = TRUE)
    return(res)
  }
  best = NMOF::DEopt(OF = perform_split,
    algo = list(min = lower, max = upper, repair = repair),
    xval = xval, y = y, min.node.size = min.node.size, objective = objective)

  return(list(split.points = sort.int(best$par), objective.value = best$value))
}

# Optimization with constrained Nelder-Mead
find_best_multiway_split_nelder = function(xval, y, n.splits = 1, min.node.size = 10,
  objective, control = NULL, ...) {
  n.splits = adjust_nsplits(xval, n.splits)
  # if only one split is needed, we use exhaustive search
  if (n.splits == 1)
    return(find_best_binary_split(xval, y, n.splits = n.splits, min.node.size = min.node.size,
      objective = objective))

  if (is.null(control))
    control = list(maxit = 640, reltol = 1e-4) #min((n.splits*2)*100, 1200)

  # find optimal split points constrained optimization
  # init = gr(rep(0, n.splits), xval, y, objective, min.node.size)
  n.quantiles = NULL #min(4*n.splits, 40)
  candidate = generate_split_candidates(xval, n.quantiles = n.quantiles, min.node.size = min.node.size)
  init = sort.int(sample(candidate, size = n.splits))
  ui = rbind(diag(n.splits), rep(0, n.splits)) - rbind(rep(0, n.splits), diag(n.splits))
  ci = c(min(xval), rep(0, n.splits - 1), -1*max(xval))
  best = constrOptim(theta = init, f = perform_split,
    method = "Nelder-Mead",
    ui = ui, ci = ci, control = control,
    xval = xval, y = y, min.node.size = min.node.size, objective = objective)

  return(list(split.points = sort.int(best$par), objective.value = best$value))
}

# Optimization with constrained simulated annealing
find_best_multiway_split_sann = function(xval, y, n.splits = 1, min.node.size = 10,
  objective, control = NULL, ...) {
  n.splits = adjust_nsplits(xval, n.splits)
  # if only one split is needed, we use exhaustive search
  if (n.splits == 1)
    return(find_best_binary_split(xval, y, n.splits = n.splits, min.node.size = min.node.size,
      objective = objective))

  # sample split points from observed x values as candidates in optimization procedure
  gr = function(i, xval, y, objective, min.node.size) {
    npar = length(i)
    q = generate_split_candidates(xval, n.quantiles = NULL, min.node.size = min.node.size)
    sort.int(sample(q, size = npar))
  }

  if (is.null(control))
    control = list(maxit = 640, reltol = 1e-4)
  # replace maximum number of evaluations depending of number of unique x values
  xval.ncomb = prod(rep(length(unique(xval)), n.splits))
  control$maxit = min(control$maxit, xval.ncomb)

  # use simulated annealing to find optimal split points
  n.quantiles = NULL #min(4*n.splits, 40)
  #init = generate_split_candidates(xval, n.quantiles = n.quantiles, min.node.size = min.node.size)[1:n.splits]
  candidate = generate_split_candidates(xval, n.quantiles = n.quantiles, min.node.size = min.node.size)
  init = sort.int(sample(candidate, size = n.splits)) #candidate[1:n.splits]
  ui = rbind(diag(n.splits), rep(0, n.splits)) - rbind(rep(0, n.splits), diag(n.splits))
  ci = c(min(xval), rep(0, n.splits - 1), -1*max(xval))
  best = constrOptim(theta = init, f = perform_split,
    method = "SANN", gr = gr,
    ui = ui, ci = ci, control = control,
    xval = xval, y = y, min.node.size = min.node.size, objective = objective)

  # init = rep(0, n.splits)
  # best = optim(init, fn = perform_split,
  #   xval = xval, y = y, min.node.size = min.node.size,
  #   objective = objective,
  #   method = "SANN", gr = gr, ...)

  return(list(split.points = sort.int(best$par), objective.value = best$value))
}
