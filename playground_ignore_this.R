nsim = 1000L
set.seed(31415)
x = sort(runif(n = nsim, min = 0, max = 5*pi))
q = quantile(x, seq(0, 1, length.out = 100), type = 1)
y = 4 + 2 * cos(x) + rnorm(nsim, mean = 0, sd = abs(cos(x)) / 2)
plot(x,y)

SS_frechet2 = function(x, y) { # slow
  # using only y-axis of curves is enough as x-axis is always the same for all curves
  require(fda.usc)
  center = unname(colMeans(y)) # this is the y-axis of pdp
  pdp.coord = cbind(center)
  dist = apply(y, 1, function(ice) Frechet(pdp.coord, cbind(ice)))
  sum(dist)
}

SS_hausdorff = function(x, y) { # slow and does not always work, buggy?
  require(SimilarityMeasures)
  center = unname(colMeans(y)) # this is the y-axis of pdp
  #grid = as.numeric(names(center)) # x-axis of pdp and ice
  #pdp.coord = cbind(grid, center)
  #dist = apply(y, 1, function(ice) metric.hausdorff(pdp.coord, cbind(grid, ice)))
  pdp.coord = cbind(center)
  dist = apply(y, 1, function(ice) metric.hausdorff(pdp.coord, cbind(ice)))
  sum(dist)
}


library(devtools)
load_all()



library(mlrMBO)


ps = makeIntegerVectorParam("split.points", len = 1L, lower = 1, upper = length(x))
#makeDiscreteVectorParam("split.points", len = 1L, values = unname(quantile(x, 1:99/100, type = 1)))
ps = makeParamSet(ps)
obj.fun = makeSingleObjectiveFunction(
  name = "objective",
  fn = perform_split,
  par.set = ps)
#ctrl = makeMBOControl()
# do three MBO iterations
#ctrl = setMBOControlTermination(ctrl, iters = 3L)
# use 500 points in the focussearch (should be sufficient for 2d)
#ctrl = setMBOControlInfill(ctrl, opt.focussearch.points = 500)
ctrl = makeMBOControl(final.method = "best.predicted", final.evals = 10)
ctrl = setMBOControlInfill(ctrl, crit = crit.ei)
ctrl = setMBOControlTermination(ctrl, iters = 7)
# create initial design
des = generateDesign(n = 5L, getParamSet(obj.fun), fun = lhs::maximinLHS)

res = mbo(obj.fun, design = des, control = ctrl, #show.info = FALSE,
  more.args	= list(xval = x, y = y, min.node.size = 1, objective = SS))

find_best_multiway_split(x, y, objective = SS)

fnMap = function(xi) sort(x[vnapply(xi, function(i) which.min(abs(x - i)))])
#.perform_split = memoise(perform_split)
n.split = 4
init = replicate(n.split, sample(x, size = n.split*10))
init = t(apply(init, 1, sort))
opt = DEoptim(fn = perform_split, lower = rep(min(x), n.split), upper = rep(max(x), n.split), x = x, y = y,
  objective = SS, min.node.size = 1, fnMap = fnMap, control = list(trace = FALSE, initialpop = init, itermax = 100))
opt$optim

opt = JDEoptim(lower = rep(min(x), n.split), upper = rep(max(x), n.split), fn = perform_split, x = x, y = y,
  objective = SS, min.node.size = 1)

o = genoud(perform_split, nvars = n.split, max = FALSE, x = x, y = y,
  objective = SS, min.node.size = 1)


.perform_split = function(par) perform_split(split.points = par, xval = x, y = y, objective = objective, min.node.size = min.node.size)
op = malschains(.perform_split, lower = rep(min(x), n.split), upper = rep(max(x), n.split), initialpop = init, verbosity = 0,
  control = malschains.control(istep = 300, ls = "sw"))
op$sol
op$fitness


# sample split points from observed x values as candidates in optimization procedure
gr = function(i, xval, y, objective, min.node.size) {
  npar = length(i)
  q = generate_split_candidates(xval, use.quantiles = use.quantiles)
  par = sample(q, size = npar)
  sort(par)
}
# use simulated annealing to find optimal split points
init = rep(0, n.split)
best = optim(init, fn = perform_split,
  xval = xval, y = y, objective = objective, min.node.size = min.node.size,
  method = "SANN", gr = gr, control = list(maxit = 10000))
best
#k = 4
#GenSA(rep(0, k), fn = perform_split, lower = rep(min(x), k), upper = rep(max(x), k), x = x, y = y, objective = SS)

library(microbenchmark)
microbenchmark(
  optim(c(0), fn = perform_split_mem, x = x, y = y, objective = SS, method = "SANN", gr = gr, control = list(trace = TRUE, maxit = 5000)),
  optim(c(0), fn = perform_split, x = x, y = y, objective = SS, method = "SANN", gr = gr, control = list(trace = TRUE, maxit = 5000)), times = 10)



# Function that finds point for numeric features
treesplit.numeric = function(x, y, objective) {
  require(BBmisc)
  # try out all points as potential split points and compute the objective in both groups
  # q = x
  q = quantile(x, seq(0, 1, length.out = 100), type = 1) # to speedup we use only quantile values as possible split points
  splits = vnapply(q, function(i) {
    left.ind = x <= i
    # split y space
    yleft = y[left.ind] # TODO: error handling if vector is empty or contains too few values
    yright = y[!left.ind] # TODO: error handling if vector is empty or contains too few values

    # split x space (not required if objective is independend from x)
    xleft = x[left.ind]
    xright = x[!left.ind]

    objective(x = xleft, y = yleft) + objective(x = xright, y = yright)
  })
  # select the split point yielding the minimal objective
  best = which.min(splits)
  split.index = which(x == q[best]) # TODO: potential bug if quantile value not observed value in x?
  list(split.index = split.index, split.point = q[best], value = splits[best])
}

# Function that finds point for numeric features
treesplit2 = function(x, y, objective) {
  require(BBmisc)
  # try out all points as potential split points and compute the objective in both groups
  # q = x
  q = quantile(x, seq(0, 1, length.out = 100), type = 1) # to speedup we use only quantile values as possible split points
  splits = vnapply(q, function(i) {
    split.factor = findInterval(x, i, rightmost.closed = TRUE) #cut2(x, i)

    y.list = split(y, split.factor)
    x.list = split(x, split.factor)

    res = vnapply(seq_along(y.list), fun = function(i) {
      objective(x = x.list[[i]], y = y.list[[i]])
    })

    sum(res)
  })
  # select the split point yielding the minimal objective
  best = which.min(splits)
  split.index = which(x == q[best]) # TODO: potential bug if quantile value not observed value in x?
  list(split.index = split.index, split.point = q[best])
}

microbenchmark(treesplit.numeric(x, y, objective = SS), treesplit2(x, y, objective = SS), times = 10)






# Determine the first split:
firstsplit_ind <- which(x <= split_lm(x,y))
lm_split1 <- lm(y[firstsplit_ind]~x[firstsplit_ind])
lm_nodeleft <- lm(y[-firstsplit_ind]~x[-firstsplit_ind])
# Visualize the split point, as well as two
# linear models on each node:
plot(x,y, cex = 0.5)
abline(v = split_lm(x,y), col = "red")
lines(x = x[firstsplit_ind],
  lm_split1$fitted.values, col = "red")
lines(x = x[-firstsplit_ind],
  lm_nodeleft$fitted.values, col = "red")

