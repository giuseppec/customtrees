nsim = 1000L
set.seed(31415)
x = sort(runif(n = nsim, min = 0, max = 5*pi))
q = quantile(x, seq(0, 1, length.out = 100), type = 1)
y = 4 + 2 * cos(x) + rnorm(nsim, mean = 0, sd = abs(cos(x)) / 2)
plot(x,y)



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




split_optimizer(x, y, objective = SS)

#opt = DEoptim(fn = perform_split, lower = c(min(x), min(x)), upper = c(max(x), max(x)), x = x, y = y, objective = SS)
#opt$optim
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

# TODO
treesplit.factor = function(x, y, objective) {
  # TODO
}

s = treesplit.numeric(x, y, objective = SS)
s

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

