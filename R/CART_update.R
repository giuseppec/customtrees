compute_variance_reductions = function(x, y, splits) {
  N = length(y)
  #mean_total = mean(y)
  #var_total = sum((y - mean_total)^2)
  
  left_indices = lapply(splits, function(split) which(x <= split))
  right_indices = lapply(splits, function(split) which(x > split))
  
  N_L = sapply(left_indices, length)
  N_R = sapply(right_indices, length)
  
  #valid_splits = N_L > 0 & N_R > 0
  
  #mean_L = sapply(left_indices, function(indices) mean(y[indices]))
  #mean_R = sapply(right_indices, function(indices) mean(y[indices]))
  
  var_L = sapply(left_indices, function(indices) sum((y[indices] - mean(y[indices]))^2))
  var_R = sapply(right_indices, function(indices) sum((y[indices] - mean(y[indices]))^2))
  
  #variance_reductions = var_total - (var_L + var_R)
  objective = var_L + var_R
  best = which.min(objective)
  
  data.frame(Split = splits[best], objective = objective[best])
  #results[valid_splits, ]
}

compute_variance_reductions1 = function(x, y, splits) {
  N = length(y)
  #mean_total = mean(y)
  #var_total = sum((y - mean_total)^2)
  
  left_indices = lapply(splits, function(split) which(x <= split))
  right_indices = lapply(splits, function(split) which(x > split))
  
  N_L = vapply(left_indices, length, FUN.VALUE = NA_integer_, USE.NAMES = FALSE)
  N_R = vapply(right_indices, length, FUN.VALUE = NA_integer_, USE.NAMES = FALSE)
  
  #valid_splits = N_L > 0 & N_R > 0
  
  #mean_L = vapply(left_indices, function(indices) mean(y[indices]), FUN.VALUE = NA_real_, USE.NAMES = FALSE)
  #mean_R = vapply(right_indices, function(indices) mean(y[indices]), FUN.VALUE = NA_real_, USE.NAMES = FALSE)
  
  var_L = vapply(left_indices, function(indices) sum((y[indices] - mean(y[indices]))^2), FUN.VALUE = NA_real_, USE.NAMES = FALSE)
  var_R = vapply(right_indices, function(indices) sum((y[indices] - mean(y[indices]))^2), FUN.VALUE = NA_real_, USE.NAMES = FALSE)

  #variance_reductions = var_total - (var_L + var_R)
  objective = var_L + var_R
  best = which.min(objective)
  
  data.frame(Split = splits[best], objective = objective[best])
  #results[valid_splits, ]
}

compute_variance_reductions2 = function(x, y, splits) {
  ord = order(x)
  x = x[ord]
  y = y[ord]
  N = length(y)
  #mean_total = mean(y)
  #var_total = sum((y - mean_total)^2)
  
  cum_sum_y = cumsum(y)
  cum_sum_y2 = cumsum(y^2)
  
  objective = vapply(splits, function(split) {
    idx = which(x <= split)
    if (length(idx) == 0 || length(idx) == N) {
      return(Inf) # Invalid split
    }
    
    N_L = length(idx)
    N_R = N - N_L
    
    S_L = cum_sum_y[N_L]
    S_R = cum_sum_y[N] - S_L
    
    SS_L = cum_sum_y2[N_L]
    SS_R = cum_sum_y2[N] - SS_L
    
    #mean_L = S_L / N_L
    #mean_R = S_R / N_R
    
    var_L = SS_L - S_L^2 / N_L
    var_R = SS_R - S_R^2 / N_R
    
    return(var_L + var_R)
  }, FUN.VALUE = NA_real_)
  
  best = which.min(objective)
  
  data.frame(Split = splits[best], objective = objective[best])
}

# sorted x
library(microbenchmark)
set.seed(1)
nsim = 5000L
x = sort(runif(n = nsim, min = 0, max = 2*pi))
q = quantile(x, seq(0, 1, length.out = 100), type = 1)
y = ifelse(x > pi/2, rnorm(nsim, mean = 0), rnorm(nsim, mean = 10, sd = 2))
psplit = generate_split_candidates(x, n.quantiles = 100, min.node.size = 10)

SS = function(y, x, requires.x = FALSE, ...) {
  ypred = mean(y)
  sum((y - ypred)^2)
}

sapply = compute_variance_reductions(x, y, psplit)
vapply = compute_variance_reductions1(x, y, psplit)
(update = compute_variance_reductions2(x, y, psplit))
(old = split_parent_node(y, data.frame(X = x), objective = SS, optimizer = find_best_binary_split))


microbenchmark(
  sapply = compute_variance_reductions(x, y, psplit), 
  vapply = compute_variance_reductions1(x, y, psplit),
  update = compute_variance_reductions2(x, y, psplit),
  #old = split_parent_node(y, data.frame(X =x), objective = SS, optimizer = find_best_binary_split),
  times = 10
  )






# unsorted x
library(microbenchmark)
set.seed(1)
nsim = 5000L
x = (runif(n = nsim, min = 0, max = 2*pi))
q = quantile(x, seq(0, 1, length.out = 100), type = 1)
y = ifelse(x > pi/2, rnorm(nsim, mean = 0), rnorm(nsim, mean = 10, sd = 2))
psplit = generate_split_candidates(x, n.quantiles = 100, min.node.size = 10)

SS = function(y, x, requires.x = FALSE, ...) {
  ypred = mean(y)
  sum((y - ypred)^2)
}

sapply = compute_variance_reductions(x, y, psplit)
vapply = compute_variance_reductions1(x, y, psplit)
(update = compute_variance_reductions2(x[order(x)], y[order(x)], psplit))
(old = split_parent_node(y, data.frame(X = x), objective = SS, optimizer = find_best_binary_split))


microbenchmark(
  sapply = compute_variance_reductions(x, y, psplit), 
  vapply = compute_variance_reductions1(x, y, psplit),
  update = compute_variance_reductions2(x, y, psplit),
  #old = split_parent_node(y, data.frame(X =x), objective = SS, optimizer = find_best_binary_split),
  times = 10
)
