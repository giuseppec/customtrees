
search_split_slow2 = function(x, Y, splits) {
  #mean_total = mean(y)
  #var_total = sum((y - mean_total)^2)
  ord = order(x)
  x = x[ord]
  Y = Y[ord,]
  N = nrow(Y)
  cum_sum_y = Rfast::colCumSums(Y)
  cum_sum_y2 = Rfast::colCumSums(Y^2)
  
  loss = rbindlist(lapply(splits, function(split) {
    idx = which(x <= split)
    if (length(idx) == 0 || length(idx) == N) {
      return(Inf) # Invalid split
    }
    
    N_L = length(idx)
    N_R = N - N_L
    
    S_L = cum_sum_y[N_L,]
    S_R = cum_sum_y[N,] - S_L
    
    SS_L = cum_sum_y2[N_L,]
    SS_R = cum_sum_y2[N,] - SS_L
    
    #mean_L = S_L / N_L
    #mean_R = S_R / N_R
    
    var_L = SS_L - S_L^2 / N_L
    var_R = SS_R - S_R^2 / N_R
    
    return(data.frame(var_L, var_R))
  }))
  
  
  loss$split.point = splits
  return(loss)
  
  
  
  risk_L = Reduce("+", lapply(risk, "[[", "var_L"))
  risk_R = Reduce("+", lapply(risk, "[[", "var_R"))
  
  objective = risk_L + risk_R
  #best = which.min(objective)
  
  data.frame(Split = splits, objective = objective, RL = risk_L, RR = risk_R)
}


search_split_slow = function(x, Y, splits) {
  ord = order(x)
  x = x[ord]
  Y = Y[ord,]
  #mean_total = mean(y)
  #var_total = sum((y - mean_total)^2)
  
  risk = lapply(Y, function(y) {
    N = length(y)
    cum_sum_y = cumsum(y)
    cum_sum_y2 = cumsum(y^2)
    
    loss = rbindlist(lapply(splits, function(split) {
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
      
      return(data.frame(var_L, var_R))
    }))
    loss$split.point = splits
    return(loss)
  })
  
  risk_L = Reduce("+", lapply(risk, "[[", "var_L"))
  risk_R = Reduce("+", lapply(risk, "[[", "var_R"))
  
  objective = risk_L + risk_R
  #best = which.min(objective)
  
  data.frame(Split = splits, objective = objective, RL = risk_L, RR = risk_R)
}
# sort something?
search_split2 = function(x, y, splits = x) {
  checkmate::assert_matrix(y)
  ord = order(x)
  x = x[ord]
  y = y[ord,]
  N = nrow(y)
  #mean_total = mean(y)
  #var_total = sum((y - mean_total)^2)
  
  cum_sum_y = Rfast::colCumSums(y)
  cum_sum_y2 = Rfast::colCumSums(y^2)
  
  loss = lapply(splits, function(split) {
    idx = which(x <= split)
    if (length(idx) == 0 || length(idx) == N) {
      return(Inf) # Invalid split
    }
    
    N_L = length(idx)
    N_R = N - N_L
    
    S_L = cum_sum_y[N_L,]
    S_R = cum_sum_y[N,] - S_L
    
    SS_L = cum_sum_y2[N_L,]
    SS_R = cum_sum_y2[N,] - SS_L
    
    #mean_L = S_L / N_L
    #mean_R = S_R / N_R
    
    var_L = (SS_L - S_L^2 / N_L)
    var_R = (SS_R - S_R^2 / N_R)
    
    return(data.frame(var_L, var_R))
  })
  
  
  risk_L = Reduce("+", lapply(loss, "[[", "var_L"))
  risk_R = Reduce("+", lapply(loss, "[[", "var_R"))
  
  objective = risk_L + risk_R
  best = which.min(objective)
  
  data.frame(Split = splits[best], objective = objective[best])
}

# 
# 
# var_per_grid = function(x, y, splits = generate_split_candidates(x, n.quantiles = 100, min.node.size = 10)) {
#   N = length(y)
#   #mean_total = mean(y)
#   #var_total = sum((y - mean_total)^2)
#   
#   cum_sum_y = cumsum(y)
#   cum_sum_y2 = cumsum(y^2)
#   
#   var = rbindlist(lapply(splits, function(split) {
#     idx = which(x <= split)
#     if (length(idx) == 0 || length(idx) == N) {
#       return(Inf) # Invalid split
#     }
#     
#     N_L = length(idx)
#     N_R = N - N_L
#     
#     S_L = cum_sum_y[N_L]
#     S_R = cum_sum_y[N] - S_L
#     
#     SS_L = cum_sum_y2[N_L]
#     SS_R = cum_sum_y2[N] - SS_L
#     
#     #mean_L = S_L / N_L
#     #mean_R = S_R / N_R
#     
#     var_L = SS_L - S_L^2 / N_L
#     var_R = SS_R - S_R^2 / N_R
#     
#     return(data.frame(var_L, var_R))
#   }))
#   var$splits = splits
#   return(var)
# }
# 
# vars = lapply(as.data.frame(Y), function(y) {
#   var_per_grid(x, y)
# })
# 
# risk_L = Reduce("+", lapply(vars, "[[", "var_L"))
# risk_R = Reduce("+", lapply(vars, "[[", "var_R"))
# 
# risk_L + risk_R
