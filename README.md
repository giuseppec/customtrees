
# Install

Install the development version from GitHub:

``` r
remotes::install_github("giuseppec/customtrees")
```

# Objectives

``` r
library(devtools)
load_all()
```

``` r
# objective that fits a constant in the nodes (CART) 
SS = function(x, y) {
  ypred = mean(y)
  sum((y - ypred)^2)
}

# objective that fits a linear model in the nodes (mob)
SS_lm = function(x, y) {
  ypred = predict(lm(y ~ x))
  sum((y - ypred)^2)
}

# objective for multivariate targets (multivariate tree), see MultivariateRandomForest::Node_cost function
SS_multi = function(x, y, cov) {
  center = colMeans(y)
  # cov = cov(y) we need to pass the cov of all data
  sum(mahalanobis(y, center = center, cov = cov, tol = 1e-30))
}
```

# Notes

  - This package is not intended to be fast. It serves as a modular
    framework and playground to explore/study the splitting of features
    by custom objectives.
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

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

# CART with multiple splits (constant model in node)

``` r
y = ifelse(x < pi/2, rnorm(nsim, mean = 0), 
  ifelse(x < pi, rnorm(nsim, mean = 10, sd = 2), 
    rnorm(nsim, mean = -10, sd = 5)))

split = split_parent_node(y, X, objective = SS, optimizer = find_best_multiway_split, 
  n.splits = 2, control = list(maxit = 5000))
split
```

    ##    feature objective.value      split.points best.split
    ## 1:       x         13897.6 1.547843,3.118802       TRUE

``` r
plot(x, y)
abline(v = unlist(split$split.points))
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

# MOB with binary splits (linear model in node)

``` r
y = 4 + 2 * cos(x) + rnorm(nsim, mean = 0, sd = abs(cos(x)) / 2)

split = split_parent_node(y, X, objective = SS_lm, optimizer = find_best_binary_split, 
  n.splits = 1, control = list(maxit = 1000))
split
```

    ##    feature objective.value split.points best.split
    ## 1:       x        140.2617     3.118802       TRUE

``` r
plot(x, y)
abline(v = unlist(split$split.points))
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

# MOB with multiple splits (linear model in node)

``` r
y = 4 + 2 * cos(x*2) + rnorm(nsim, mean = 0, sd = abs(cos(x)) / 2)

split = split_parent_node(y, X, objective = SS_lm, optimizer = find_best_multiway_split, 
  n.splits = 3, control = list(maxit = 1000))
split
```

    ##    feature objective.value               split.points best.split
    ## 1:       x        171.9338 1.480383,3.118802,4.529151       TRUE

``` r
plot(x, y)
abline(v = unlist(split$split.points))
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

# Group ICE Curves with Multivariate Tree (binary splits, constant model in node)

We first generate some functional data:

``` r
# Simulate Data
n = 1000
x1 = runif(n, -1, 1)
x2 = runif(n, -1, 1)
x3 = sample(c(0, 1), size = n, replace = TRUE, prob = c(0.5, 0.5))
x4 = sample(c(0, 1), size = n, replace = TRUE, prob = c(0.7, 0.3))
eps = rnorm(n, 0, 1)

# noisy vars
x5 = sample(c(0, 1), size = n, replace = TRUE, prob = c(0.5, 0.5))
x6 = rnorm(n, mean = 1, sd = 5)

y = 0.2*x1 - 8*x2 + ifelse(x3 == 0, I(16*x2),0) + ifelse(x1 > mean(x1), I(8*x2),0) + eps

dat = data.frame(x1, x2, x3, x4, x5, x6, y)
X = dat[, setdiff(colnames(dat), "y")]

# Fit model and compute ICE for x2
library(iml)
library(ranger)
mod = ranger(y ~ ., data = dat, num.trees = 10)
pred = predict.function = function(model, newdata) predict(model, newdata)$predictions
model = Predictor$new(mod, data = X, y = dat$y, predict.function = pred)
effect = FeatureEffects$new(model, method = "ice", grid.size = 20, features = "x2")

# Plot ICE curves: WE WANT TO FIND SUBGROUPS SUCH THAT ICE KURVES ARE HOMOGENOUS
library(ggplot2)
ggplot(effect$results$x2, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id))
```

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

Formulate curves above by multivariate target and find feature that
splits the curves such that they are more homogenous in the nodes:

``` r
# Get ICE values and arrange them in a horizontal matrix
library(tidyverse)
Y = spread(effect$results$x2, .borders, .value)
Y = Y[, setdiff(colnames(Y), c(".type", ".id", ".feature"))]
str(X) # contains our feature values
```

    ## 'data.frame':    1000 obs. of  6 variables:
    ##  $ x1: num  -0.155 0.322 -0.849 0.637 0.789 ...
    ##  $ x2: num  -0.862 -0.624 -0.188 0.782 -0.969 ...
    ##  $ x3: num  0 1 1 0 1 1 1 0 1 1 ...
    ##  $ x4: num  0 1 0 0 0 1 0 0 0 0 ...
    ##  $ x5: num  0 1 0 0 0 1 0 1 1 0 ...
    ##  $ x6: num  5.24 7.98 1.01 -3.7 -5.82 ...

``` r
str(Y) # contains ICE values for each grid point
```

    ## 'data.frame':    1000 obs. of  20 variables:
    ##  $ -0.996347735170275 : num  -8.19 0.285 4.818 -10.952 0.666 ...
    ##  $ -0.891456903546656 : num  -7.401 0.126 4.753 -10.323 0.487 ...
    ##  $ -0.786566071923038 : num  -6.4942 0.0211 5.3642 -10.3731 0.2382 ...
    ##  $ -0.681675240299419 : num  -5.116 -0.243 4.865 -9.421 -1.015 ...
    ##  $ -0.576784408675801 : num  -4.438 -0.943 3.133 -7.499 -1.434 ...
    ##  $ -0.471893577052182 : num  -4.015 -0.943 3.482 -5.781 -0.21 ...
    ##  $ -0.367002745428564 : num  -3.909 -0.122 1.448 -5.193 -0.22 ...
    ##  $ -0.262111913804945 : num  -3.383 -0.433 1.407 -4.561 -0.347 ...
    ##  $ -0.157221082181327 : num  -2.081 -0.284 0.928 -1.485 0.456 ...
    ##  $ -0.0523302505577081: num  -1.734 -0.284 1.006 -0.564 0.297 ...
    ##  $ 0.0525605810659104 : num  -1.0066 -0.0526 -0.0878 0.9909 0.025 ...
    ##  $ 0.157451412689529  : num  0.0866 -0.019 -0.8669 1.7487 -0.4287 ...
    ##  $ 0.262342244313147  : num  2.056 -0.609 -1.779 3.401 -1.162 ...
    ##  $ 0.367233075936766  : num  2.221 -0.178 -1.955 5.374 -0.499 ...
    ##  $ 0.472123907560384  : num  2.221 -0.157 -2.615 5.946 -0.365 ...
    ##  $ 0.577014739184003  : num  2.799 -0.358 -3.281 6.553 -0.389 ...
    ##  $ 0.681905570807622  : num  4.236 0.865 -3.185 8.982 -0.313 ...
    ##  $ 0.78679640243124   : num  4.236 0.695 -5.354 9.306 -1.331 ...
    ##  $ 0.891687234054859  : num  4.803 0.374 -4.644 10.949 -0.321 ...
    ##  $ 0.996578065678477  : num  4.79 2.27 -3.63 10.24 1.34 ...

``` r
# compute covariance for data and use this in for mahalanobis distance in the objective
COV = cov(Y)
SS_multi2 = function(x, y) SS_multi(x, y, cov = COV)
sp = split_parent_node(Y = Y, X = X, objective = SS_multi2, 
  n.splits = 1, min.node.size = 1, optimizer = find_best_binary_split)
sp
```

    ##    feature objective.value split.points best.split
    ## 1:      x1        19373.03 -0.008479368      FALSE
    ## 2:      x2        19937.16   0.05241069      FALSE
    ## 3:      x3        19066.63          0.5       TRUE
    ## 4:      x4        19282.92          0.5      FALSE
    ## 5:      x5        19475.43          0.5      FALSE
    ## 6:      x6        19505.03     9.783206      FALSE

``` r
node_index = generate_node_index(Y, X, result = sp)
str(node_index)
```

    ## List of 2
    ##  $ 0.0      : int [1:480] 1 4 8 13 14 16 17 20 21 22 ...
    ##  $ [0.5,1.0]: int [1:520] 2 3 5 6 7 9 10 11 12 15 ...

``` r
# Compare with MultivariateRandomForest yields same result
library(MultivariateRandomForest)
invcov = solve(cov(Y), tol = 1e-30)
sp2 = splitt2(X = as.matrix(X), Y = as.matrix(Y), m_feature = ncol(X), 
  Index = 1:nrow(X), Inv_Cov_Y = invcov, Command = 2, ff = 1:ncol(X))
str(sp2)
```

    ## List of 4
    ##  $ Idx_left       : int [1:480] 1 4 8 13 14 16 17 20 21 22 ...
    ##  $ Idx_right      : int [1:520] 2 3 5 6 7 9 10 11 12 15 ...
    ##  $ Feature_number : int 3
    ##  $ Threshold_value: num 0.5

Visualize the results:

``` r
plot.data = effect$results$x2
plot.data$.split = as.factor(ifelse(plot.data$.id %in% node_index[[1]], "left.node", "right.node"))

ggplot(plot.data, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id)) + facet_grid(~ .split)
```

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->
