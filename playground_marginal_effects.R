
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
library(partykit)
library(ggparty)
load_all()


tsk = makeRegrTask(data = Boston, target = "medv")
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
  model = mod, feature = "rm", data = Boston[, !colnames(Boston) %in% c("marginal.effect", "medv")], step.size = 2, objective = SS,
  target.ratio.sd.mean = 0.4, min.node.size = 30)

leaf.nodes.list = lapply(return.list, FUN = function(x) {
  if (is.na(x["split.feature"])) {
    return(x[c("id", "data")])
  } else {
    return(NULL)
  }
})

leaf.nodes.list = leaf.nodes.list[!unlist(lapply(leaf.nodes.list, is.null))]

tree.summary = tree_summary(return.list)
leaf.nodes = get_leaf_nodes(tree.summary)
party.tree = party_interface(tree.summary, Boston)
split.sequence.largest = record_split_sequence(tree.summary, "10")
feature.subspace.largest = get_feature_subspace(Boston, split.sequence.largest)
get_global_subspace(Boston, leaf.nodes)

p = ggparty(party.tree) +
  geom_edge() +
  geom_edge_label() +
  geom_node_splitvar() +
  geom_node_info() +
  geom_node_plot(gglist = list(geom_boxplot(aes(x = marginal.effect)),
                               xlab("Marginal Effects"),
                               theme_bw(),
                               coord_flip(),
                               theme(axis.text.x = element_blank(),
                                     axis.ticks.x = element_blank(),
                                     axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)))),
                               nudge_y = -0.11,
                               shared_axis_labels = TRUE)
                 
p

leaf.nodes
subspace = get_global_subspace(data = Boston, leaf.nodes = leaf.nodes)
# leaf.largest = subspace[as.numeric(as.character(subspace$n.observations)) == max(as.numeric(as.character(subspace$n.observations))), ]
# 
# save(leaf.largest, file = "ame_subspace_largest.Rdata")
save(subspace, file = "ame_subspace.Rdata")
ggsave(file = "ggparty.png", plot = p, width = 10, height = 8, dpi = 150, units = "in", device='png')
