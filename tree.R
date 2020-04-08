nsim = 1000L
set.seed(31415)
x = sort(runif(n = nsim, min = 0, max = 5*pi))
q = quantile(x, seq(0, 1, length.out = 100), type = 1)
y = 4 + 2 * cos(x) + rnorm(nsim, mean = 0, sd = abs(cos(x)) / 2)
plot(x,y)

# objective for CART
SS = function(x, y) {
  ypred = mean(y)
  sum((y - ypred)^2)
}

# objective for mob
SS_mod <- function(x, y){
  ypred = predict(lm(y ~ x))
  sum((y - ypred)^2)
}

# Split function
# treesplit = function(x, y, objective) {
#   UseMethod("treesplit")
# }

# input: x feature values, y target, min.node.size, n.split.points
optFun = function(split.point, x, y, min.node.size = 10, objective = SS) {
  split.factor = findInterval(x, split.point, rightmost.closed = TRUE)
  node.size = table(split.factor)

  if (any(node.size < min.node.size))
    return(Inf)

  d = data.table(x, y, split.factor)
  sum(d[, .(obj = objective(x, y)), by = split.factor]$obj)
}


optFun2 = function(split.point, x, y, min.node.size = 10, objective = SS_mod) {
  split.point = sort(split.point)
  split.factor = findInterval(x, split.point, rightmost.closed = TRUE)
  node.size = table(split.factor)

  if (any(node.size < min.node.size))
    return(Inf)

  y.list = split(y, split.factor)
  x.list = split(x, split.factor)

  res = vnapply(1:length(y.list), fun = function(i) {
    objective(x = x.list[[i]], y = y.list[[i]])
  })

  sum(res)
}
GenSA(rep(0, 3), fn = optFun2, lower = rep(min(x), 3), upper = rep(max(x), 3), x = x, y = y)

k = 4
gr = function(i, x, y) sort(sample(quantile(x, seq(0, 1, length.out = 100), type = 1), size = k))
a = optim(rep(0, k), fn = optFun2, x = x, y = y, method = "SANN", gr = gr, control = list(trace = TRUE))
a
microbenchmark(optFun(q[5], x, y, objective = SS), optFun2(q[5], x, y, objective = SS), times = 100)



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

    res = vnapply(1:length(y.list), fun = function(i) {
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

split_lm <- function(x, y) {
  # try out all points as potential split points and ...
  splits <- lapply(x, function(i) {
    yleft <- y[x <= i]
    yright <- y[x > i]
    #SS(yleft) + SS(yright) # ... compute SS in both groups
    # Create a subset of values of x that correspond
    # to the values of y, i.e. ---> (x,y) \in N_j
    xleft <- x[x <= i]
    xright <- x[x > i]
    if (length(xright) < 1) {
      SS_mod(resp = yleft, feature = xleft)
    } else {
      SS_mod(resp = yleft, feature = xleft) + SS_mod(resp = yright, feature = xright)
    }
  })
  # select the split point yielding the minimal sum of squares (SS)
  best <- which.min(splits)
  x[best]
}
split_lm(x, y) # the first split point:




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

