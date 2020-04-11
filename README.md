
# Install

To make use of functions from this package, you need to clone this
repository, install the `devtools` R package, navigate to the directory
of this package and use the `load_all()` function.

``` r
# install.packages("devtools")
library(devtools)
load_all()
```

Install the development version from GitHub:

``` r
# Currently, documentation is missing. Hence, installation may fail
# remotes::install_github("giuseppec/customtrees")
```

# Notes

  - This package is not intended to be fast. It serves as a modular
    framework and playground to explore/study the splitting of features
    by custom objectives.
  - Currently only trees of depth 1 are fitted with `split_parent_node`.
    If you want a tree, you need to call this function recursively on
    the generated child nodes. You can use `generate_node_index` to get
    the split indices of the observations from the current node.
  - Splits for categorical variables currently not implemented and
    tested. Try to handle categoricals as numerics as workaround.
  - The `perform_split` function computes (and aggregates) the objective
    in the generated nodes after splitting w.r.t. specific split points.
  - Binary splits generate two nodes and are implemented in
    `find_best_binary_split`. The implementation does exhaustive search
    of split point candidates to find the best split point for a given
    feature.
  - Multiple splits generate multiple nodes and are implemented in
    `find_best_multiway_split`. The implementation currently uses a slow
    simulated annealing optimization to find the best split point for a
    given feature (might be improved and replaced with other, faster
    optimization procedures).

# Define Objectives used as Split Criteria

``` r
library(tidyverse)
library(Rmalschains)
library(iml)
library(ranger)
library(kmlShape)
library(dtw)
```

``` r
# objective that fits a constant in the nodes (CART) 
SS = function(y, x, requires.x = FALSE) {
  ypred = mean(y)
  sum((y - ypred)^2)
}

# objective that fits a linear model in the nodes (mob)
SS_lm = function(y, x, requires.x = TRUE) {
  ypred = predict(lm(y ~ x))
  sum((y - ypred)^2)
}

# objective for multivariate targets (multivariate tree), see MultivariateRandomForest::Node_cost function
SS_mah = function(y, x, requires.x = FALSE, cov) {
  center = colMeans(y)
  # cov = cov(y) we need to pass the cov of all data
  sum(mahalanobis(y, center = center, cov = cov, tol = 1e-30))
}

# Frechet distance FDA measure
SS_fre = function(y, x, requires.x = FALSE) { # slow
  # using only y-axis of curves is enough as x-axis is always the same for all curves
  require(kmlShape)
  center = colMeans(y)
  grid.x = as.numeric(names(center))
  pdp.y = unname(center)
  dist = apply(y, 1, function(ice) distFrechet(grid.x, pdp.y, grid.x, ice, FrechetSumOrMax = "sum"))
  sum(dist)
}

# Dynamic time warping FDA measure
SS_dtw = function(y, x, requires.x = FALSE) {
  require(dtw)
  pdp = colMeans(y) # this is the pdp
  dist = apply(y, 1, function(ice) dtw(ice, pdp, distance.only = TRUE)$normalizedDistance)
  sum(dist)
}
```

# CART with binary splits (constant model in node)

``` r
nsim = 1000L
x = x = sort(runif(n = nsim, min = 0, max = 2*pi))
q = quantile(x, seq(0, 1, length.out = 100), type = 1)
y = ifelse(x > pi/2, rnorm(nsim, mean = 0), rnorm(nsim, mean = 10, sd = 2))
X = data.frame(x = x)

split = split_parent_node(y, X, objective = SS, optimizer = find_best_binary_split)
split
```

    ##    feature objective.value split.points best.split
    ## 1:       x        1900.012     1.547843       TRUE

``` r
# plot result
plot(x, y)
abline(v = unlist(split$split.points))
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

# Extending CART to multiple splits (constant model in node)

``` r
y = ifelse(x < pi/2, rnorm(nsim, mean = 0), 
  ifelse(x < pi, rnorm(nsim, mean = 10, sd = 2), 
    rnorm(nsim, mean = -10, sd = 5)))

# Simulated annealing from find_best_multiway_split
split = split_parent_node(y, X, objective = SS, optimizer = find_best_multiway_split, 
  n.splits = 2, control = list(maxit = 1000))

# MA-LS Chains optimization from find_best_multiway_split2 (claimed to be faster and better)
split2 = split_parent_node(y, X, objective = SS, optimizer = find_best_multiway_split2, 
  n.splits = 2, control = malschains.control(istep = 100, ls = "sw"))

split
```

    ##    feature objective.value      split.points best.split
    ## 1:       x        12695.16 1.571378,3.142373       TRUE

``` r
split2
```

    ##    feature objective.value      split.points best.split
    ## 1:       x        12695.16 1.570233,3.140220       TRUE

``` r
plot(x, y)
abline(v = unlist(split$split.points))
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

# MOB with binary splits (linear model in node)

``` r
y = 4 + 2 * cos(x) + rnorm(nsim, mean = 0, sd = abs(cos(x)) / 2)

split = split_parent_node(y, X, objective = SS_lm, optimizer = find_best_binary_split, n.splits = 1)
split
```

    ##    feature objective.value split.points best.split
    ## 1:       x        148.9909     3.118802       TRUE

``` r
plot(x, y)
abline(v = unlist(split$split.points))
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

# MOB with multiple splits (linear model in node)

``` r
y = 4 + 2 * cos(x*2) + rnorm(nsim, mean = 0, sd = abs(cos(x)) / 2)

split = split_parent_node(y, X, objective = SS_lm, optimizer = find_best_multiway_split2, 
  n.splits = 3, control = malschains.control(istep = 100, ls = "sw"))
split
```

    ##    feature objective.value               split.points best.split
    ## 1:       x         150.879 1.574027,3.166505,4.711258       TRUE

``` r
plot(x, y)
abline(v = unlist(split$split.points))
```

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

# Group ICE Curves with Multivariate Tree (binary splits, constant model in node)

We first generate some functional data:

``` r
# Simulate Data
n = 500
x1 = runif(n, -1, 1)
x2 = runif(n, -1, 1)
x3 = sample(c(0, 1), size = n, replace = TRUE, prob = c(0.5, 0.5))
x4 = sample(c(0, 1), size = n, replace = TRUE, prob = c(0.7, 0.3))
eps = rnorm(n, 0, 1)

# noisy vars
x5 = sample(c(0, 1), size = n, replace = TRUE, prob = c(0.5, 0.5))
x6 = rnorm(n, mean = 1, sd = 5)

y = 0.2*x1 - 8*x2 + ifelse(x3 == 0, I(16*x2),0) + ifelse(x1 > mean(x1), I(8*x2),0) + eps
# We also get interesting results using a 2-way interaction of numeric features
# y = 0.2*x1 - 8*x2 + 8*x6*x2 + eps
# y = 0.2*x1 - 8*x2^2 + 5*cos(x2*5)*x6 + ifelse(x3 == 0, I(8*x2),0) + eps

dat = data.frame(x1, x2, x3, x4, x5, x6, y)
X = dat[, setdiff(colnames(dat), "y")]

# Fit model and compute ICE for x2
mod = ranger(y ~ ., data = dat, num.trees = 500)
pred = predict.function = function(model, newdata) predict(model, newdata)$predictions
model = Predictor$new(mod, data = X, y = dat$y, predict.function = pred)
effect = FeatureEffects$new(model, method = "ice", grid.size = 20, features = "x2")

# Plot ICE curves: WE WANT TO FIND SUBGROUPS SUCH THAT ICE KURVES ARE HOMOGENOUS
ggplot(effect$results$x2, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id))
```

![](README_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

Formulate curves above by multivariate target and find feature that
splits the curves such that they are more homogenous in the nodes:

``` r
# Get ICE values and arrange them in a horizontal matrix
Y = spread(effect$results$x2, .borders, .value)
Y = Y[, setdiff(colnames(Y), c(".type", ".id", ".feature"))]
str(X) # contains our feature values
```

    ## 'data.frame':    500 obs. of  6 variables:
    ##  $ x1: num  -0.647 0.265 -0.577 -0.184 0.395 ...
    ##  $ x2: num  -0.321 -0.428 0.615 -0.731 -0.29 ...
    ##  $ x3: num  1 1 1 0 1 1 0 0 0 1 ...
    ##  $ x4: num  0 0 0 0 0 1 1 0 0 0 ...
    ##  $ x5: num  1 0 1 0 1 0 1 1 1 1 ...
    ##  $ x6: num  -4.967 0.803 0.387 7.317 -3.323 ...

``` r
str(Y) # contains ICE values for each grid point
```

    ## 'data.frame':    500 obs. of  20 variables:
    ##  $ -0.999021098483354 : num  4.598 0.626 4.363 -5.492 1.293 ...
    ##  $ -0.894629949140117 : num  3.623 0.42 3.462 -5.754 0.424 ...
    ##  $ -0.790238799796881 : num  3.714 0.373 3.663 -5.002 0.353 ...
    ##  $ -0.685847650453644 : num  3.5862 0.0551 3.5012 -5.3176 0.1411 ...
    ##  $ -0.581456501110408 : num  3.2422 -0.2229 3.094 -5.0661 -0.0586 ...
    ##  $ -0.477065351767171 : num  2.799 -0.324 2.309 -4.036 0.199 ...
    ##  $ -0.372674202423935 : num  2.506 -0.314 1.826 -3.25 0.298 ...
    ##  $ -0.268283053080698 : num  1.767 -0.182 1.29 -2.492 0.402 ...
    ##  $ -0.163891903737462 : num  1.1459 -0.0796 0.5469 -1.5278 0.3551 ...
    ##  $ -0.0595007543942254: num  0.406 -0.154 -0.25 -1.026 0.248 ...
    ##  $ 0.0448903949490111 : num  -0.169 0.101 -0.929 0.184 0.466 ...
    ##  $ 0.149281544292248  : num  -0.9643 -0.0779 -1.6312 0.8887 0.1218 ...
    ##  $ 0.253672693635484  : num  -1.6837 -0.0956 -2.1791 2.4092 -0.0932 ...
    ##  $ 0.358063842978721  : num  -2.02248 -0.04181 -2.41665 3.17457 0.00383 ...
    ##  $ 0.462454992321957  : num  -2.3342 0.0613 -2.8936 3.8548 0.0787 ...
    ##  $ 0.566846141665194  : num  -2.64 -0.0973 -3.2978 4.1874 -0.0584 ...
    ##  $ 0.67123729100843   : num  -3.1794 -0.1359 -3.6819 4.464 -0.0227 ...
    ##  $ 0.775628440351667  : num  -3.889 0.314 -4.086 4.966 0.198 ...
    ##  $ 0.880019589694903  : num  -4.5693 -0.0503 -4.6865 4.7836 -0.0414 ...
    ##  $ 0.98441073903814   : num  -3.986 -0.438 -4.055 4.51 0.174 ...

``` r
# compute covariance for data and use this in for mahalanobis distance in the objective
COV = cov(Y)
SS_mah2 = function(y, x, requires.x = FALSE) 
  SS_mah(y = y, x = x, requires.x = requires.x, cov = COV)
sp = split_parent_node(Y = Y, X = X, objective = SS_mah2, 
  n.splits = 1, optimizer = find_best_binary_split)
sp
```

    ##    feature objective.value split.points best.split
    ## 1:      x1        9621.242 -0.009225438      FALSE
    ## 2:      x2        9931.309   -0.2088459      FALSE
    ## 3:      x3        9502.743          0.5       TRUE
    ## 4:      x4        9599.929          0.5      FALSE
    ## 5:      x5        9570.492          0.5      FALSE
    ## 6:      x6        9605.409     1.818777      FALSE

``` r
node_index = generate_node_index(Y, X, result = sp)
str(node_index)
```

    ## List of 2
    ##  $ class: Factor w/ 2 levels "[0,0.5]","(0.5,1]": 2 2 2 1 2 2 1 1 1 2 ...
    ##  $ index:List of 2
    ##   ..$ [0,0.5]: int [1:253] 4 7 8 9 11 13 14 15 16 17 ...
    ##   ..$ (0.5,1]: int [1:247] 1 2 3 5 6 10 12 19 24 28 ...

``` r
# frechet distance yields same splits but is a bit slower
sp_frechet = split_parent_node(Y = Y, X = X, objective = SS_fre,
  n.splits = 1, optimizer = find_best_binary_split)
sp_frechet
```

    ##    feature objective.value split.points best.split
    ## 1:      x1       25581.779   -0.8686654      FALSE
    ## 2:      x2       25816.928  -0.00974361      FALSE
    ## 3:      x3        9049.959          0.5       TRUE
    ## 4:      x4       25932.858          0.5      FALSE
    ## 5:      x5       25998.220          0.5      FALSE
    ## 6:      x6       25765.512     10.03504      FALSE

``` r
node_index_frechet = generate_node_index(Y, X, result = sp_frechet)
str(node_index_frechet)
```

    ## List of 2
    ##  $ class: Factor w/ 2 levels "[0,0.5]","(0.5,1]": 2 2 2 1 2 2 1 1 1 2 ...
    ##  $ index:List of 2
    ##   ..$ [0,0.5]: int [1:253] 4 7 8 9 11 13 14 15 16 17 ...
    ##   ..$ (0.5,1]: int [1:247] 1 2 3 5 6 10 12 19 24 28 ...

``` r
## dynamic time warping distance yields same splits but is much slower
# sp_dtw = split_parent_node(Y = Y, X = X, objective = SS_dtw, 
#   n.splits = 1, optimizer = find_best_binary_split)
# sp_dtw
# node_index_dtw = generate_node_index(Y, X, result = sp_dtw)
# str(node_index_dtw)
```

``` r
# Compare with MultivariateRandomForest yields same result
library(MultivariateRandomForest)
invcov = solve(cov(Y), tol = 1e-30)
sp2 = splitt2(X = as.matrix(X), Y = as.matrix(Y), m_feature = ncol(X), 
  Index = 1:nrow(X), Inv_Cov_Y = invcov, Command = 2, ff = 1:ncol(X))
str(sp2)
```

    ## List of 4
    ##  $ Idx_left       : int [1:253] 4 7 8 9 11 13 14 15 16 17 ...
    ##  $ Idx_right      : int [1:247] 1 2 3 5 6 10 12 19 24 28 ...
    ##  $ Feature_number : int 3
    ##  $ Threshold_value: num 0.5

Visualize the results:

``` r
plot.data = effect$results$x2
plot.data$.split = node_index$class[plot.data$.id]

ggplot(plot.data, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id)) + facet_grid(~ .split)
```

![](README_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

# Group ICE Curves with Multivariate Tree (multiway splits, constant model in node)

Multiway split **fails** with `SS_mah2` (mahalanobis distance) as
objective. This is because the curve structure along the x-axis is not
considered in the distance calculation\!

``` r
sp_multiway = split_parent_node(Y = Y, X = X, objective = SS_mah2, 
  n.splits = 2, optimizer = find_best_multiway_split2)
sp_multiway
```

    ##    feature objective.value                split.points best.split
    ## 1:      x1        9337.197 -0.4445809963, 0.0006376415       TRUE
    ## 2:      x2        9888.553       -0.2282386, 0.3992555      FALSE
    ## 3:      x3        9502.743                         0.5      FALSE
    ## 4:      x4        9599.929                         0.5      FALSE
    ## 5:      x5        9570.492                         0.5      FALSE
    ## 6:      x6        9349.277       -0.4254728, 4.8104691      FALSE

``` r
node_index_multiway = generate_node_index(Y, X, result = sp_multiway)
str(node_index_multiway)
```

    ## List of 2
    ##  $ class: Factor w/ 3 levels "[-0.996,-0.445]",..: 1 3 1 2 3 3 3 3 3 3 ...
    ##  $ index:List of 3
    ##   ..$ [-0.996,-0.445]  : int [1:138] 1 3 12 14 24 30 34 35 36 37 ...
    ##   ..$ (-0.445,0.000638]: int [1:125] 4 11 13 16 17 18 19 21 23 31 ...
    ##   ..$ (0.000638,0.994] : int [1:237] 2 5 6 7 8 9 10 15 20 22 ...

``` r
plot.data$.split = node_index_multiway$class[plot.data$.id]

ggplot(plot.data, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id)) + facet_grid(~ .split)
```

![](README_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

Instead, using a distance measure that is suited for curves (e.g.,
frechet distance) works:

``` r
sp_multiway_frechet = split_parent_node(Y = Y, X = X, objective = SS_fre, 
  n.splits = 2, optimizer = find_best_multiway_split2)
sp_multiway_frechet
```

    ##    feature objective.value            split.points best.split
    ## 1:      x1       25105.463   -0.9498307,-0.8806692      FALSE
    ## 2:      x2       25541.096 0.001550442,0.026824208      FALSE
    ## 3:      x3        9049.959                     0.5       TRUE
    ## 4:      x4       25932.858                     0.5      FALSE
    ## 5:      x5       25998.220                     0.5      FALSE
    ## 6:      x6       25492.118       10.06241,11.09778      FALSE

``` r
node_index_multiway_frechet = generate_node_index(Y, X, result = sp_multiway_frechet)
str(node_index_multiway_frechet)
```

    ## List of 2
    ##  $ class: Factor w/ 2 levels "[0,0.5]","(0.5,1]": 2 2 2 1 2 2 1 1 1 2 ...
    ##  $ index:List of 2
    ##   ..$ [0,0.5]: int [1:253] 4 7 8 9 11 13 14 15 16 17 ...
    ##   ..$ (0.5,1]: int [1:247] 1 2 3 5 6 10 12 19 24 28 ...

``` r
plot.data$.split = node_index_multiway_frechet$class[plot.data$.id]

ggplot(plot.data, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id)) + facet_grid(~ .split)
```

![](README_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

## Group ICE Curves with Multivariate Tree (multiway splits, constant model in node)

Now, we try a non-linear effect with a continous interaction effect.

``` r
y = 0.2*x1 - 8*x2^2 + 5*cos(x2*5)*x6 + ifelse(x3 == 0, I(8*x2),0) + eps

dat = data.frame(x1, x2, x3, x4, x5, x6, y)
X = dat[, setdiff(colnames(dat), "y")]

# Fit model and compute ICE for x2
mod = ranger(y ~ ., data = dat, num.trees = 1000)
pred = predict.function = function(model, newdata) predict(model, newdata)$predictions
model = Predictor$new(mod, data = X, y = dat$y, predict.function = pred)
effect = FeatureEffects$new(model, method = "ice", grid.size = 20, features = "x2")

# Plot ICE curves: WE WANT TO FIND SUBGROUPS SUCH THAT ICE KURVES ARE HOMOGENOUS
ggplot(effect$results$x2, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id))

sp_multi = lapply(1:5, function(i) {
  split_parent_node(Y = Y, X = X, objective = SS_fre, 
  n.splits = i, optimizer = find_best_multiway_split2, min.node.size = 10)
})

rbindlist(sp_multi, idcol = "n.splits")

node_index_multiway_frechet = generate_node_index(Y, X, result = sp_multi[[5]])
str(node_index_multiway_frechet)

plot.data$.split = node_index_multiway_frechet$class[plot.data$.id]

ggplot(plot.data, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id)) + facet_grid(~ .split)
```
