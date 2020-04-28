

get_leaf_nodes = function(tree.summary) {
  
  leaf.nodes = tree.summary[is.na(tree.summary$child.left) & is.na(tree.summary$child.right), ]
  return(leaf.nodes) 
}

tree_summary = function(tree) {
  sub.list = lapply(return.list, FUN = function(x) {
    x[c("id", "depth", "parent.node", "split.feature", "split.value", "child.left", "child.right", "ame", "ame.sd", "n.observations")]
  })

  tree.summary = data.frame(do.call(rbind, sub.list))
  tree.summary = data.frame(apply(tree.summary, MARGIN = 2, FUN = unlist))
 
  tree.summary$depth = as.numeric(as.character(tree.summary$depth))
  tree.summary$split.value = as.numeric(as.character(tree.summary$split.value))
  tree.summary$ame = round(as.numeric(as.character(tree.summary$ame)), 2)
  tree.summary$ame.sd = round(as.numeric(as.character(tree.summary$ame.sd)), 2)
  
  
  return(tree.summary)
}
