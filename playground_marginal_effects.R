
library(mlr)
library(devtools)
library(tidyverse)
library(Rmalschains)
library(iml)
library(ranger)
library(kmlShape)
library(dtw)
library(mlr)
library(printr)
library(MASS)
library(dplyr)
library(tidyr)
library(tidyverse)
load_all()

data = Boston
select(data)
tsk = makeRegrTask(data = data, target = "medv")
mod = train("regr.randomForest", tsk)

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

objective.function = function(x, y, requires.x = FALSE) {
  sd(y)
}

marginal_effect_tree(model = mod, feature = "rm", data = data,
                                step.size = 1, SS)
return.list

child.nodes = lapply(return.list, FUN = function(x) {
  if (is.null(x[["child.left"]]) && is.null(x[["child.right"]])) {
    return(x)}
})
length()
node.sizes = lapply(child.nodes, FUN = function(x) {
  x[["n.observations"]]
})
summary(unlist(node.sizes))

child.nodes.cumulative.sd = sum(unlist(lapply(child.nodes, FUN = function(x) x[["ame.sd"]])))

sum(unlist(lapply(child.nodes, FUN = function(x) x[["n.observations"]]))) == nrow(Boston)

child.nodes

unlist(lapply(child.nodes, FUN = function(x) x["ame"]))

sub.list = lapply(return.list, FUN = function(x) {
  x[c("id", "split.feature", "split.value", "child.left", "child.right", "ame", "ame.sd", "n.observations")]
})

tree.summary = as.data.frame(do.call(rbind, sub.list))

