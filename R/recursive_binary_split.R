recursive_binary_split_marginal_effects = function(
  model, feature, step.size, parent.node, return.list, objective, parent.index) {
  
  marginals = marginal_effects(model, parent.node$data, feature, step.size)
  split.result = split_parent_node(
    Y = marginals, X = parent.node$data, objective = objective,
    optimizer = find_best_binary_split)
  
  child.nodes = create_child_nodes(parent.node, split.result)
  child.node.1 = child.nodes[["child.1"]]
  child.node.2 = child.nodes[["child.2"]]
  
  parent.node = set_node_split.feature(
    parent.node, get_split_data(split.result)[["split.feature"]])
  parent.node = set_node_split.value(
    parent.node, get_split_data(split.result)[["split.value"]])
  parent.node = set_node_child.left(parent.node, child.node.1$id)
  parent.node = set_node_child.right(parent.node, child.node.2$id)
  parent.node = set_node_ame(parent.node, mean(marginals))
  
  return.list = append(return.list, parent.node)
  
  marginals.child.1 = marginal_effects(
    model, child.node.1$data, feature, step.size)
  marginals.child.2 =  marginal_effects(
    model, child.node.2$data, feature, step.size)
  
  child.node.1 = set_node_ame(child.node.1, mean(marginals.child.1))
  child.node.2 = set_node_ame(child.node.1, mean(marginals.child.2))
  
  if (sd(marginals.child.1) >= 0.5 * abs(mean(marginals.child.1)) && length(marginals.child.1) > 2) {
    return(recursive_binary_split_marginal_effects(
      model = model, feature = feature, step.size = step.size,
      parent.node = child.node.1, return.list = return.list,
      objective = objective))
  } else {
    return.list = append(return.list, child.node.1)
  }
  
  if (sd(marginals.child.2) >= 0.5 * abs(mean(marginals.child.1)) && length(marginals.child.2) > 2) {
    return(recursive_binary_split_marginal_effects(
      model = model, feature = feature, step.size = step.size,
      parent.node = child.node.2, return.list = return.list,
      objective = objective))
  } else {
    return.list = append(return.list, child.node.2)
  }
  return(return.list)
}
