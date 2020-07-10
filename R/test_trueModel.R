# test with true model

source("R/perform_split.r")
source("R/split_optimizer.r")
source("R/split_parent_node.r")
source("R/helper_functions.r")
source("R/detect_interactions.r")

# 1. Extrapolation Example

# Simulate Data
set.seed(1234)

z <- runif(500,0,1)
d = ifelse(z > 0.5,1,0)
eps = rnorm(500, 0, 0.15)

# noisy vars
x5 = sample(c(0, 1), size = 500, replace = TRUE, prob = c(0.5, 0.5))
x6 = rnorm(500, mean = 1, sd = 5)


y = z^0.5 * (1 - d) + (2^(1.5*z) - 1) * d + eps
dat = data.frame(z,d,x5,x6,y)
X = dat[, setdiff(colnames(dat), "y")]
#write.csv2(dat, "data/ExtrapolationExample1.csv")

#return(result)
library(BBmisc)
library(MASS)

mymodel = makeS3Obj("mymodel", fun = function(data) return(data$z^0.5 * (1 - data$d) + (2^(1.5*data$z) - 1) * data$d))
predict.mymodel = function(object, newdata) {
  p = object$fun(newdata)
  p + rnorm(length(p), mean = 0, sd = 0.05)
}

test = predict(mymodel, dat)
mod = Predictor$new(model = mymodel, data = X, predict.function = predict.mymodel)
effect = FeatureEffects$new(predictor = mod, method = "ice", grid.size = 20, features = "z")
effect$plot()
effectd = FeatureEffects$new(predictor = mod, method = "ice", grid.size = 20, features = "x5")
effectd$plot()

# calculate interactions
potential.interactions = find_potential_interactions(X, mod, find_best_multiway_split2, SS_fre_filtered, SS_fre, 1, min.node.size = 30, improve.first.split = 0, improve.n.splits = 0)
interactions = potential.interactions[order(potential.interactions$improvement, decreasing = TRUE),]

# compare to HStatistics
ia1 = Interaction$new(mod, feature = "z", grid.size = 20)
ia2 = Interaction$new(mod, feature = "d", grid.size = 20)
ia3 = Interaction$new(mod, feature = "x5", grid.size = 20)
ia4 = Interaction$new(mod, feature = "x6", grid.size = 20)

#-------------------------------------------------------------------------------------------------------------
# simulation studies with correlated features

source("R/helper_functions_sim.r")

# generate data
p = 6
n = 500

set.seed(1234)
V1 = rnorm(n)
V2 = rnorm(n)
V3 = rnorm(n)
V4 = rnorm(n)
V5 = rnorm(n)
V6 = rnorm(n)
Xsim = data.frame(V1,V2,V3,V4,V5,V6)

Xsim = sim_data(N = n, n_groups = 1, y = 1, size_groups = 6, prop = c(0.5))
epsilon = rnorm(n)

Xsim = as.data.frame(Xsim)
colnames(Xsim) = paste0("V",1:p)

#y = Xsim[,2] + 5*Xsim[,3]^3*Xsim[,4] + ifelse(Xsim[,6] <= mean(Xsim[,6]), I(8*Xsim[,3]),0) + epsilon
y = 0.2*Xsim[,2] + 5*cos(5*Xsim[,3])*Xsim[,4] + ifelse(Xsim[,6] <= mean(Xsim[,6]), I(8*Xsim[,3]),0) + epsilon


dat = data.frame(Xsim,y)
X = dat[, setdiff(colnames(dat), "y")]

#mymodel = makeS3Obj("mymodel", fun = function(data) return(data$V2 + 5*data$V3^3*data$V4 + ifelse(data$V6 <= mean(data$V6), I(8*data$V3),0)))
mymodel = makeS3Obj("mymodel", fun = function(data) return(0.2*data$V2 + 5*cos(5*data$V3)*data$V4 + ifelse(data$V6 <= mean(data$V6), I(8*data$V3),0)))
predict.mymodel = function(object, newdata) {
  p = object$fun(newdata)
  p + rnorm(length(p), mean = 0, sd = 0.1)
}


# interactions # todo: Fehler bei "ind" nochmal testen
potential.interactions = find_potential_interactions(X, mod, find_best_multiway_split2, SS_fre_filtered, SS_fre, 1, min.node.size = 30, improve.first.split = 0, improve.n.splits = 0)
interactions = potential.interactions[order(potential.interactions$improvement, decreasing = TRUE),]
#find_true_interactions(X, mod, find_best_multiway_split2, n.splits = 1)



#----------------------------------------------------------------------------------
# Plot splitted ice curves 

test = predict(mymodel, dat)
mod = Predictor$new(model = mymodel, data = X, predict.function = predict.mymodel)
effect = FeatureEffects$new(predictor = mod, method = "ice", grid.size = 20, features = "V4")
effect$plot()

# Get ICE values and arrange them in a horizontal matrix
Y = spread(effect$results$V4, .borders, .value)
Y = Y[, setdiff(colnames(Y), c(".type", ".id", ".feature"))]
# centered ICE curves
for (i in 1:nrow(Y)) {
  Y[i,] = as.numeric(unname(Y[i,])) - mean(as.numeric(unname(Y[i,])))
}


# Plot ICE curves
ggplot(effect$results$V4, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id))


# splitted ice curves
res = split_parent_node(Y = Y, X = X, feat = "V4", optimizer = find_best_multiway_split2, n.splits = 5, objective = SS_fre_filtered, min.node.size = 10)
ice = get_ice_curves(X = X, Y = Y, result = res, extrapol = FALSE)


plot.data = plot.prep.full(ice)

# Plot splitted ICE curves
ggplot(plot.data, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id)) + facet_grid(~ .split)







