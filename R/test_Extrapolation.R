
# simulation example to test  functions
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


# Fit model and compute ICE for z
mod = ranger(y ~ ., data = dat, num.trees = 500)
pred = predict.function = function(model, newdata) predict(model, newdata)$predictions
model = Predictor$new(mod, data = X, y = y, predict.function = pred)
effect = FeatureEffects$new(model, method = "ice", grid.size = 20, features = "z")
# Get ICE values and arrange them in a horizontal matrix
Y = spread(effect$results$z, .borders, .value)
Y = Y[, setdiff(colnames(Y), c(".type", ".id", ".feature"))]
# centered ICE curves
for (i in 1:nrow(Y)) {
  Y[i,] = as.numeric(unname(Y[i,])) - mean(as.numeric(unname(Y[i,])))
}


# Plot ICE curves
ggplot(effect$results$z, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id))


#------------------------------------------------------------------------------------------------------------
# test plotting with extrapolation

res = split_parent_node(Y = Y, X = X, feat = "z", optimizer = find_best_multiway_split2, n.splits = 1, objective = SS_fre_filtered, min.node.size = 10)
ice = get_ice_curves(X = X, Y = Y, result = res)


plot.data = plot.prep.full(ice)

# Plot splitted ICE curves
ggplot(plot.data, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id)) + facet_grid(~ .split)



#------------------------------------------------------------------------------------------------------------
# test interaction functions


potential.interactions = find_potential_interactions(X, model, find_best_multiway_split2, SS_fre_filtered, SS_fre, 3, min.node.size = 10, improve.first.split = 0.1, improve.n.splits = 0.1)
interactions = find_true_interactions(X, model, find_best_multiway_split2)

# compare with hstatistics
ia1 = Interaction$new(model, feature = "z", grid.size = 20)
ia2 = Interaction$new(model, feature = "d", grid.size = 20)
ia3 = Interaction$new(model, feature = "x5", grid.size = 20)
ia4 = Interaction$new(model, feature = "x6", grid.size = 20)
ia = Interaction$new(model, grid.size = 20)




#-------------------------------------------------------------------------------------------------------------
# simulation studies with correlated features

source("R/helper_functions_sim.r")

# generate data
p = 6
n = 500

set.seed(1234)
Xsim = sim_data(N = n, n_groups = 1, y = 1, size_groups = 6, prop = c(1))
epsilon = rnorm(n)

Xsim = as.data.frame(Xsim)
colnames(Xsim) = paste0("V",1:p)

y = Xsim[,2] + 5*cos(Xsim[,3]*5)*Xsim[,4] + ifelse(Xsim[,6] <= mean(Xsim[,6]), I(8*Xsim[,3]),0) + epsilon


dat = data.frame(Xsim,y)
X = dat[, setdiff(colnames(dat), "y")]



# Fit model and compute ICE for z
# task = makeRegrTask(data = dat, target = "y")
# res = resample("regr.ranger", task, cv5)
mod = ranger(y ~ ., data = dat, num.trees = 1000)
pred = predict.function = function(model, newdata) predict(model, newdata)$predictions
model = Predictor$new(mod, data = X, y = y, predict.function = pred)

# interactions # todo: Fehler bei "ind" nochmal testen
potential.interactions = find_potential_interactions(X, model, find_best_multiway_split2, SS_fre_filtered, SS_fre, 3, min.node.size = 30, improve.first.split = 0.1, improve.n.splits = 0.1)
find_true_interactions(X, model, find_best_multiway_split2)
interactions = potential.interactions[order(potential.interactions$improvement, decreasing = TRUE),]
#saveRDS(interactions, "interactions.rds")
# zu hohe improvements durchweg - gefilterte objective funktioniert hier nicht - kann ungefilterte einfach 
# verwendet werden? nochmals durchdenken
# Interaktion zwischen x3 und x10 wird eindeutig gefunden. Nicht eindeutig für x3 und x6 (auch nicht bei HStatistik)
# Evtl. mit angepasster objective

# compare with hstatistics
ia = Interaction$new(model, feature = colnames(X)[2], grid.size = 20)




effect = FeatureEffects$new(model, method = "ice", grid.size = 20, features = "V3")
# Get ICE values and arrange them in a horizontal matrix
Y = spread(effect$results$V3, .borders, .value)
Y = Y[, setdiff(colnames(Y), c(".type", ".id", ".feature"))]
# centered ICE curves
for (i in 1:nrow(Y)) {
  Y[i,] = as.numeric(unname(Y[i,])) - mean(as.numeric(unname(Y[i,])))
}


res = split_parent_node(Y, X, feat = "V3", n.splits = 1, min.node.size = 30, optimizer = find_best_multiway_split2, extrapol = TRUE, objective = SS_fre_filtered)

ice = get_ice_curves(X = X, Y = Y, result = res)


plot.data = plot.prep.full(ice)

# Plot splitted ICE curves
ggplot(plot.data, aes(x = .borders, y = .value)) + 
  geom_line(aes(group = .id)) + facet_grid(~ .split)
