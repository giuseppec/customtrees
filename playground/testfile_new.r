library(devtools)
load_all()

source("playground/helper_functions.r")


# test differences between L1 and L2

# objective with L2 loss (with L1 loss is defined in helper_functions.r, sourced above)
SS_fre_L2 = function(y, x, requires.x = FALSE, ...) { # slow
  #require(Rfast)
  ypred = apply(y, 2, mean) #Rfast::colMedians(as.matrix(y))
  sum(t((t(y) - ypred)^2))
}


# generate data  
p = 6
n = 100

set.seed(1234)
for (i in 1:100) {
  
  V1 = rnorm(n)
  V2 = rnorm(n)
  V3 = rnorm(n)
  V4 = rnorm(n, 1, 3)
  V5 = rnorm(n)
  V6 = rnorm(n)
  Xsim = data.frame(V1,V2,V3,V4,V5,V6)
  #V5[sample(1:n, 15)] = rnorm(15, 20,10)
  
  epsilon = rnorm(n)
  #y = 0.2*Xsim[,2] + 5*Xsim[,3]*Xsim[,4] + ifelse(Xsim[,6] <= mean(Xsim[,6]), I(8*Xsim[,3]),0) + epsilon
  y = 0.2*Xsim[,2] + ifelse(Xsim[,4] <= mean(Xsim[,4]), I(5*Xsim[,3]),0) + ifelse(Xsim[,6] <= mean(Xsim[,6]), I(8*Xsim[,3]),0) + epsilon
  
  dat = data.frame(Xsim,y)
  X = dat[, setdiff(colnames(dat), "y")]
  
  
  # define model and predictor function
  mod = ranger(y ~ ., data = dat, num.trees = 1000)
  pred = predict.function = function(model, newdata) predict(model, newdata)$predictions
  model = Predictor$new(mod, data = X, y = y, predict.function = pred)
  
  
  
  # extract ICE curves - here e.g. for X3
  effect = FeatureEffect$new(model, method = "ice", grid.size = 20, feature = "V3")
  # Get ICE values and arrange them in a horizontal matrix
  colnames(effect$results)[1] = ".borders"
  Y = spread(effect$results, .borders, .value)
  Y = Y[, setdiff(colnames(Y), c(".type", ".id"))]
  for (j in 1:nrow(Y)) {
    Y[j,] = as.numeric(unname(Y[j,])) - mean(as.numeric(unname(Y[j,])))
  }
  
  # artificially create randomly for a couple of ice curves outliers (play here a bit around to show the
  # difference between L1 and L2 loss)
  Y[sample(1:100, 6),] = c(7*cos(rnorm(20, 1, 4)) - 3*as.numeric(colnames(Y)), 7*cos(rnorm(20, 2, 2)) - (as.numeric(colnames(Y)))^2,
                           5*cos(rnorm(20, 3, 2)) - 10*as.numeric(colnames(Y)), 10*sin(rnorm(20, 0, 2)) - 6*as.numeric(colnames(Y)),
                           #6*cos(rnorm(20, -1, 2)) - 4*as.numeric(colnames(Y)), 4*sin(rnorm(20, -1, 3)) - 2*as.numeric(colnames(Y)),
                           3*(rnorm(20, -2, 2))^2 - as.numeric(colnames(Y)), (rnorm(20, 0, 2))^3 - 2*as.numeric(colnames(Y)))

  # here: binary split, you can also try more splits, but I think for the showcase binary would be enough
  # here objective currently set to SS_fre which is equivalent to L1 loss, 
  # change it to SS_fre_L2 to calculate the same with L2 loss
  res.L1 = split_parent_node(Y = Y, X = X[,-which(colnames(X) == "V3")], 
                          optimizer = find_best_multiway_split, n.splits = 1, objective = SS_fre, min.node.size = 10,
                          feat = "V3", x.all = X)
  res.L1$rel.improve = 1 - (res.L1$objective.value/SS_fre(Y,X))
  res.L1 = res.L1[order(res.L1$objective.value),]
  res.L1$rank = 1:(p - 1)
  res.L2 = split_parent_node(Y = Y, X = X[,-which(colnames(X) == "V3")], 
                             optimizer = find_best_multiway_split, n.splits = 1, objective = SS_fre_L2, min.node.size = 10,
                             feat = "V3", x.all = X)
  res.L2$rel.improve = 1 - (res.L2$objective.value/SS_fre_L2(Y,X))
  res.L2 = res.L2[order(res.L2$objective.value),]
  res.L2$rank = 1:(p - 1)
  
  if (i == 1) {
    result.total.L1 = res.L1
    result.total.L2 = res.L2
  }
  else {
    result.total.L1 = rbind(result.total.L1, res.L1)
    result.total.L2 = rbind(result.total.L2, res.L2)
  }
   
}

# overview on Ranking and relative improvement
aggregate(result.total.L1[,6:7], list(result.total.L1$feature), summary)
aggregate(result.total.L2[,6:7], list(result.total.L2$feature), summary)



# overview on split points for L1 - for interaction features V6 and V4
summary(unlist(result.total.L1$split.points[which(result.total.L1$feature == "V6")]))
summary(unlist(result.total.L1$split.points[which(result.total.L1$feature == "V4")]))


# overview on split points for L2 - for interaction features V6 and V4
summary(unlist(result.total.L2$split.points[which(result.total.L1$feature == "V6")]))
summary(unlist(result.total.L2$split.points[which(result.total.L1$feature == "V4")]))




# plot ice curves
ice = get_ice_curves(X = X, Y = Y, result = result.total.L1, extrapol = FALSE)
plot.data = plot.prep.full(ice)

ggplot(plot.data, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id)) + facet_grid(~ .split)

# play a bit around, compare it visually and the output of "res" --> try to find explanations: 
# which cases give similar results and which don't
# if you think you found a good example to compare L1 and L2, then do
# 100 iterations and extract splitting features and splitting points (maybe also ranks)
# L1 loss should be more robust / consistent


