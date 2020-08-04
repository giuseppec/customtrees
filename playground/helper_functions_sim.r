library(mlr)
make_correlated_feats = function(n, rv, name, prop) {
  for (i in 1:n) {
    if (i == 1) {
      U_i = rv
      random_20 = sample(1:length(U_i), prop*length(U_i))
      #U_i[random_20] = rnorm(10, 0, 1)
      U_i[random_20] = 0.5*rv[random_20] + 0.5*rnorm(prop*length(U_i), 0, 1)
      df = setNames(data.frame(U_i), paste0(name, i))
    } else {
      U_i = rv
      random_20 = sample(1:length(U_i), prop*length(U_i))
      U_i[random_20] = 0.5*rv[random_20] + 0.5*rnorm(prop*length(U_i), 0, 1)
      df = data.frame(df, setNames(data.frame(U_i), paste0(name, i)))
    }
  }
  df
}

sim_data = function(N = 100, n_groups, y = numeric(), size_groups = numeric(), prop = numeric()) {
  x = rnorm(N, 0, 1)
  for (i in 1:length(size_groups)) {
    
    if (i == 1) {
      df = make_correlated_feats(size_groups[i], rv = x, name = paste0("Group_", i, "_feat_"), prop = prop[i])
      rv = list(x)
      next
    }
    df = data.frame(df, make_correlated_feats(size_groups[i], rv = x, name = paste0("Group_", i, "_feat_"), prop = prop[i]))
    rv = c(rv, list(x))
  }
  checkmate::assert_true(length(y) == n_groups)
  for (i in 1:length(y)) {
    if (i == 1) {
      a = list(y[i] * rv[[i]])
      z = a[[1]]
      next
    }
    a = c(a, list(y[i] * rv[[i]]))
    z = z + a[[i]]
  }
  #target =  setNames(data.frame(z - mean(z) + rnorm(N, 0, 0.1)), "target")
  data.frame(df)
}
