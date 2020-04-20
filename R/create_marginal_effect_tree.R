create_marginal_effect_tree = function(model, feature, data, step.size, objective) {
  
  assign("return.list", list(), envir = .GlobalEnv)
  start.node = Node(id = "0", data = data)
  iter.count = 0
  assign("iter.count", iter.count, envir = .GlobalEnv)
  
  recursive_binary_split_marginal_effects(
    model = model, feature = feature, step.size = step.size,
    this.node = start.node, objective = objective)
}

get_leaf_nodes = function(tree.summary) {
  
  leaf.nodes = tree.summary[is.na(tree.summary$child.left) & is.na(tree.summary$child.right), ]
  return(leaf.nodes) 
}

tree_summary = function(tree) {
  sub.list = lapply(return.list, FUN = function(x) {
    x[c("id", "split.feature", "split.value", "child.left", "child.right", "ame", "ame.sd", "n.observations")]
  })

  tree.summary = data.frame(do.call(rbind, sub.list))
  tree.summary = data.frame(apply(tree.summary, MARGIN = 2, FUN = unlist))
  tree.summary$id = as.character(tree.summary$id)
  return(tree.summary)
}
