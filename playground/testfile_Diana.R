library(devtools)
library(BBmisc)
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
V1 = rnorm(n)
V2 = rnorm(n)
V3 = rnorm(n, 2, 6)
V4 = rnorm(n)
V5 = rnorm(n, mean = 3, sd = 5)
V6 = rnorm(n, mean = 3, sd = 6)
Xsim = data.frame(V1,V2,V3,V4,V5,V6)

epsilon = rnorm(n)
y = 0.2*Xsim[,2] + 5*Xsim[,3]*Xsim[,4] + ifelse(Xsim[,6] <= mean(Xsim[,6]), I(8*Xsim[,3]),0) + epsilon

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
Y[sample(1:100, 4),] = c(10*cos(rnorm(20, 1, 5)) - 3*as.numeric(colnames(Y)), 7*cos(rnorm(20, 2, 4)) - (as.numeric(colnames(Y)))^2,
  5*cos(rnorm(20, 3, 3)) - 20*as.numeric(colnames(Y)), 20*sin(rnorm(20, 0, 2)) - 8*as.numeric(colnames(Y)))

# here: binary split, you can also try more splits, but I think for the showcase binary would be enough
# here objective currently set to SS_fre which is equivalent to L1 loss, 
# change it to SS_fre_L2 to calculate the same with L2 loss
res = split_parent_node(Y = Y, X = X[,-which(colnames(X) == "V3")], 
  optimizer = find_best_multiway_split, n.splits = 1, objective = SS_fre, min.node.size = 10,
  feat = "V3", x.all = X)

# plot ice curves
ice = get_ice_curves(X = X, Y = Y, result = res, extrapol = FALSE)
plot.data = plot.prep.full(ice)

ggplot(plot.data, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id)) + facet_grid(~ .split)

# play a bit around, compare it visually and the output of "res" --> try to find explanations: 
# which cases give similar results and which don't
# if you think you found a good example to compare L1 and L2, then do
# 100 iterations and extract splitting features and splitting points (maybe also ranks)
# L1 loss should be more robust / consistent


#-------------------------------------------------------------------------------

# interaction detection functions

# for the ranking simulation you can also use the following functions:
# they will calculate all two-sided interactions between all features which
# fulfill that the best splitting feature first split improvement is at least 10% 
# (improvement of binary split compared to root node).
# number of splits depends on the relative improvement which has to be at least
# improve.n.splits. The number is limited by n.splits. You can play here around as well

rm(res)
potential.interactions = find_potential_interactions(x = X, model = model, optimizer = find_best_multiway_split, objective.split = SS_fre, objective.total = SS_fre, 
  n.splits = 5, min.node.size = 10, improve.first.split = 0.1, improve.n.splits = 0.1, 
  x.all = X)
interactions = potential.interactions[order(potential.interactions$improvement, decreasing = TRUE),]



# for the ranking example please try using the real underlying model as discussed before
# adjust the formula in "mymodel" according to your simulation example
# you can play around with the sd of the randomness you add to the parameters 
# (otherwise you have no variation in your ice curves)
# you can use the "model" as you did before with the random forest

mymodel = makeS3Obj("mymodel", fun = function(data) return(data$V2 + 5*data$V3^3*data$V4 + ifelse(data$V6 <= mean(data$V6), I(8*data$V3),0)))
predict.mymodel = function(object, newdata) {
  p = object$fun(newdata)
  p + rnorm(length(p), mean = 0, sd = 0.1)
}
model = Predictor$new(model = mymodel, data = X, predict.function = predict.mymodel)



