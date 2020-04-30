
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
mod = train("regr.ksvm", tsk)

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

SS_variation_coeff = function(y, x, requires.x = FALSE) {
  sd(y) / abs(mean(y))
}

SS_quartile_disp_coeff = function(y, x, requires.x = FALSE) {
  diff(quantile(y, c(0.25, 0.75))) / median(y)
}

objective.function = function(x, y, requires.x = FALSE) {
  sd(y)
}

marginal.effects = marginal_effects(mod, Boston, "rm", 1)
Boston$marginal.effect = marginal.effects


set.seed(123)

marginal_effect_tree(
  model = mod, feature = "rm",
  data = Boston[, !colnames(Boston) %in% c("marginal.effect", "medv")],
  step.size = 1, objective = SS_variation_coeff,
  target.objective = 0.5, min.node.size = 50)

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

party.tree[[1]]$data$marginal.effect

tree.summary$variation.coefficient = (tree.summary$ame.sd / abs(tree.summary$ame))
p = ggparty(party.tree, add_vars = list("ame" =
                                          function(data, node) {
                                            mean(node$data$marginal.effect)
                                          },
                                        "sd" =
                                          function(data, node) {
                                            sd(node$data$marginal.effect)
                                          },
                                        "variation.coefficient" =
                                          function(data, node) {
                                          sd(node$data$marginal.effect) / abs(mean(node$data$marginal.effect))
                                        },
                                        "target.objective" =
                                          function(data, node) {
                                            sd(node$data$marginal.effect) / abs(mean(node$data$marginal.effect)) < 0.5
                                          })) +
  geom_edge() +
  geom_edge_label() +
  geom_node_splitvar() +
  # geom_node_info() + 
  geom_node_label(aes(label = info, color = target.objective), ids = "terminal", size = 3, nudge_y = -0.1) +
  geom_node_plot(gglist = list(geom_boxplot(aes(x = marginal.effect)),
                               xlab("Marginal Effects"),
                               theme_bw(),
                               coord_flip(),
                               theme(axis.text.x = element_blank(),
                                     axis.ticks.x = element_blank())),
                 nudge_y = -0.15,
                 shared_axis_labels = TRUE) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) +
  scale_color_discrete("Objective target value achieved") +
  ggtitle("Subspace Average Marginal Effects of Feature 'rm' with Step Size 1")

p


subspace = get_global_subspace(data = Boston, leaf.nodes = leaf.nodes)
# leaf.largest = subspace[as.numeric(as.character(subspace$n.observations)) == max(as.numeric(as.character(subspace$n.observations))), ]
# 
# save(leaf.largest, file = "ame_subspace_largest.Rdata")
save(subspace, file = "ame_subspace.Rdata")
ggsave(file = "ggparty.png", plot = p, width = 16, height = 8, dpi = 300, units = "in", device='png')
