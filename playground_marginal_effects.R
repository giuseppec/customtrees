
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

tree.summary = tree_summary(return.list)
leaf.nodes = get_leaf_nodes(tree.summary)

# child.nodes = lapply(return.list, FUN = function(x) {
  # if (is.null(x[["child.left"]]) && is.null(x[["child.right"]])) {
    # return(x)}
# })
# length()
# node.sizes = lapply(child.nodes, FUN = function(x) {
  # x[["n.observations"]]
# })
# summary(unlist(node.sizes))
# 
# child.nodes.cumulative.sd = sum(unlist(lapply(child.nodes, FUN = function(x) x[["ame.sd"]])))
# 
# sum(unlist(lapply(child.nodes, FUN = function(x) x[["n.observations"]]))) == nrow(Boston)
# 
# child.nodes
# 
# unlist(lapply(child.nodes, FUN = function(x) x["ame"]))
# 
# sub.list = lapply(return.list, FUN = function(x) {
#   x[c("id", "split.feature", "split.value", "child.left", "child.right", "ame", "ame.sd", "n.observations")]
# })

# tree.summary = data.frame(do.call(rbind, sub.list))
# tree.summary = data.frame(apply(tree.summary, MARGIN = 2, FUN = unlist))
# # tree.summary$id = as.character(tree.summary$id)
# 
# tree.leaf.nodes = tree.summary[is.na(tree.summary$split.feature), ]

# feature.subspaces = lapply(leaf.nodes$id, FUN = function(id) {
#   split.sequence = record_split_sequence(tree = tree.summary, leaf.id = id)
#   feature.subspace.list = get_feature_subspace(data, split.sequence)
#   return(feature.subspace.list)
# })
# 
# feature.subspaces = lapply(feature.subspaces, FUN = function(feature.subspace) {
#   complementary.features = colnames(data)[!colnames(data) %in% colnames(feature.subspace)]
#   feature.subspace[complementary.features] = NA 
#   return(feature.subspace)
# })

    
leaf.nodes
subspace = get_global_subspace(data = data, leaf.nodes = leaf.nodes)
head(subspace)
save(subspace, file = "ame_subspace.Rdata")
