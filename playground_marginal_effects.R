
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

marginal.effects = marginal_effects(mod, Boston, "rm", 1)
Boston$marginal.effect = marginal.effects

set.seed(123)
marginal_effect_tree(
  model = mod, feature = "rm", data = Boston[, !colnames(Boston) %in% "marginal.effect"], step.size = 2, objective = SS,
  target.ratio.sd.mean = 0.4)
return.list
leaf.nodes.list = lapply(return.list, FUN = function(x) {
  if (is.na(x["split.feature"])) {
    return(x[c("id", "data")]) 
  } else {
    return(NULL)
  }
})

leaf.nodes.list = leaf.nodes.list[!unlist(lapply(leaf.nodes.list, is.null))]
# drop NULL values
tree.summary = tree_summary(return.list)
leaf.nodes = get_leaf_nodes(tree.summary)

####
party_interface = function(tree.summary, data) {
  options(scipen = 99)
  
  features = colnames(data)
  party.node.list = list()
  depth = max(tree.summary$depth)
  nodes = tree.summary[which(tree.summary$depth == depth), ]
  nodes.id = as.numeric(as.character(nodes$id))
  
  party.nodes = lapply(nodes.id, FUN = function(node.id) {
    
    node.ame = nodes[nodes$id == node.id, "ame"]
    node.sd = nodes[nodes$id == node.id, "ame.sd"]
    node.n = nodes[nodes$id == node.id, "n.observations"]
    node.frac = round(as.numeric(as.character(node.n)) / nrow(data) * 100, 2)
    node.frac = paste0(node.frac, "%")
    node.info = paste0("AME: ", node.ame, "\n", "SD: ", node.sd, "\n", "Frac.: ", node.frac)
    party.node = partynode(id = node.id, info = node.info)
    return(party.node)
  })
  
  party.node.list = append(party.node.list, party.nodes)
  depth = depth - 1
  
  while (depth > 0) {
    
    nodes = tree.summary[tree.summary$depth == depth, ]
    parent.nodes = nodes[!is.na(nodes$split.feature), ]
    leaf.nodes = nodes[is.na(nodes$split.feature), ]
    
    parent.party.nodes = lapply(1:nrow(parent.nodes), FUN = function(i) {
      node = parent.nodes[i, ]
      node.id = as.numeric(as.character(node$id))
      node.split.feature = node$split.feature
      node.split.value = node$split.value
      node.split.feature.index = which(features == node.split.feature)
      
      child.party.nodes = lapply(party.node.list, FUN = function(x) {
        
        x.splitstring = unlist(strsplit(as.character(x$id), ""))
        x.parent = paste0(x.splitstring[-length(x.splitstring)], collapse = "")
        if (node.id == x.parent) {
          return(x)
        } else {
        }
      })
      
      child.party.nodes = child.party.nodes[!unlist(lapply(child.party.nodes, is.null))]
      # drop NULL values
      
      party.node = partynode(
        node.id, split = partysplit(node.split.feature.index, node.split.value),
        kids = child.party.nodes)
      return(party.node)
    })
    
    if (nrow(leaf.nodes) > 0) {
      leaf.party.nodes = lapply(1:nrow(leaf.nodes), FUN = function(i) {
        node = leaf.nodes[i, ]
        node.id = as.numeric(as.character(node$id))
        node.ame = nodes[nodes$id == node.id, "ame"]
        node.sd = nodes[nodes$id == node.id, "ame.sd"]
        node.n = nodes[nodes$id == node.id, "n.observations"]
        node.frac = round(as.numeric(as.character(node.n)) / nrow(data) * 100, 2)
        node.frac = paste0(node.frac, "%")
        node.info = paste0("AME: ", node.ame, "\n", "SD: ", node.sd, "\n", "Frac.: ", node.frac)
        party.node = partynode(id = node.id, info = node.info)
        return(party.node)
      })
    } else {}
    
    party.node.list = append(party.node.list, parent.party.nodes)
    party.node.list = append(party.node.list, leaf.party.nodes)
    
    depth = depth - 1
  }
  
  tree.index = which(unlist(lapply(party.node.list, FUN = function(x) {
    x$id == 1
  })))
  
  party.tree = party.node.list[[tree.index]]
  party.tree = party(party.tree, data)
                     # fitted = data.frame(
                     #   "(fitted)" = fitted_node(party.tree,
                     #                            data = data),
                     #   "(response)" = data$marginal.effect),
                     # terms = terms(marginal.effect ~ ., data = data),)
  return(party.tree)
}

leaf.nodes
party.tree = party_interface(tree.summary, Boston)
plot(party.tree)  
# p = autoplot(party.tree) +
#   geom_edge() +
#   geom_edge_label() +
#   geom_node_splitvar()

party.tree[10]
ggplot(party.tree[10]$data) +
  geom_boxplot(aes(x = marginal.effect))

p = ggparty(party.tree) +
  geom_edge() +
  geom_edge_label() +
  geom_node_splitvar() +
  # geom_node_info(fill = "green") + 
  geom_node_plot(gglist = list(geom_boxplot(aes(x = marginal.effect))))
p

ggplot(data = party.tree[21]$data) +
  geom_boxplot(aes(x = ""))

ggplot(data = Boston[1:2,]) +
  geom_boxplot(aes(x = marginal.effect))

# 
# nc <- length(unique(ggplot_build(p)[["data"]][[1]][["PANEL"]]))
# 
# ggsave("tall.pdf", p, width=7, height = nc * 20 + 2, limitsize = FALSE)


theme_update(text = element_text(size=0.5))
p
# p + theme(axis.text.x = element_text(color = "grey20", size = 1, angle = 0, hjust = .5, vjust = .5, face = "plain"),
#         axis.text.y = element_text(color = "grey20", size = 1, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
#         axis.title.x = element_text(color = "grey20", size = 1, angle = 0, hjust = .5, vjust = 0, face = "plain"),
#         axis.title.y = element_text(color = "grey20", size = 1, angle = 0, hjust = .5, vjust = .5, face = "plain"))
  
# node = tree.summary[tree.summary$id == nodes.parent, ]
  # node.id = as.numeric(as.numeric(as.character(node$id)))
  # depth = node$depth
  
  node.split.feature = node$split.feature
  node.split.feature.id = which(features == node.split.feature)
  node.split.value = node$split.value
  
  party.node.1 = partynode(node.id, split = partysplit(node.split.feature.id, breaks = node.split.value),
            kids = party.nodes)
  
  #
  
  
  ####
  
  
  

####
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
leaf.largest = subspace[as.numeric(as.character(subspace$n.observations)) == max(as.numeric(as.character(subspace$n.observations))), ]
# 
# save(leaf.largest, file = "ame_subspace_largest.Rdata")
# save(subspace, file = "ame_subspace.Rdata")
