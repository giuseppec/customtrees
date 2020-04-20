
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

tree.summary = data.frame(do.call(rbind, sub.list))
tree.summary = data.frame(apply(tree.summary, MARGIN = 2, FUN = unlist))
tree.summary$id = as.character(tree.summary$id)

tree.leaf.nodes = tree.summary[is.na(tree.summary$split.feature), ]

record_split_sequence = function(tree, leaf.id) {
  
  id.split = unlist(strsplit(x = leaf.id, split = "", fixed = TRUE))
  split.list = list()
  id = character(1)
  
  for (i in 1:(length(id.split) - 1)) {
    
    id = paste0(id, id.split[i])
    split = tree[tree$id == id, c("split.feature", "split.value")]
    if (id.split[i + 1] == "0") {
      split$direction = "left" 
    } else {
      split$direction = "right"  
    }
    split.list = append(split.list, list(split))
  }
  split.sequence = data.frame(do.call(rbind, split.list))
  split.sequence$split.feature = as.character(split.sequence$split.feature)
  split.sequence$split.value = as.numeric(as.character(split.sequence$split.value))
  return(split.sequence)
}
get_feature_subspace = function(data, split.sequence) {
  
  split.features = unique(split.sequence$split.feature)
  
  feature.subspace.list = lapply(split.features, FUN = function(feature) {
    feature.splits = split.sequence[split.sequence$split.feature == feature, ]
    
    if ("left" %in% feature.splits$direction) {
      boundary.right = min(feature.splits$split.value[feature.splits$direction == "left"])
    } else {
      boundary.right = max(data[[feature]])
    }
    if ("right" %in% feature.splits$direction) {
      boundary.left = max(feature.splits$split.value[feature.splits$direction == "right"])
    } else {
      boundary.left = min(data[[feature]])
    }
    feature.subspace = paste0("Feature '", feature, "': [", boundary.left, ", ", boundary.right, "]")
    return(feature.subspace)
  })
  return(feature.subspace.list)
}

# split.sequence = record_split_sequence(tree = tree.summary, leaf.id = leaf.id)
# feature.subspace.list = get_feature_subspace(data, split.sequence)
# 

tree.leaf.nodes

lapply(tree.leaf.nodes$id, FUN = function(id) {
  split.sequence = record_split_sequence(tree = tree.summary, leaf.id = id)
  feature.subspace.list = get_feature_subspace(data, split.sequence)
  return(feature.subspace.list)
})

