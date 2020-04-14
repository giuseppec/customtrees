recursive_binary_split_marginal_effects = function(
  model, feature, step.size, parent.node, return.list, objective) {
  
  marginals = marginal_effects(model, parent.node, feature, step.size)
  split.result = split_parent_node(
    Y = marginals, X = parent.node, objective = objective,
    optimizer = find_best_binary_split)
  
  if (nrow(split.result) >= 1) {
    child.nodes = create_child_nodes(parent.node, split.result)
    child.node.1 = child.nodes[[1]]
    child.node.2 = child.nodes[[2]]
    marginals.child.1 = marginal_effects(
      model, child.node.1, feature, step.size)
    print(length(marginals.child.1))
    marginals.child.2 =  marginal_effects(
      model, child.node.2, feature, step.size)
    print(length(marginals.child.2))
    
    if (sd(marginals.child.1) >= 0.5 * abs(mean(marginals.child.1)) && length(marginals.child.1) > 2) {
      return.list = append(recursive_binary_split_marginal_effects(
        model = model, feature = feature, step.size = step.size,
        parent.node = child.node.1, return.list = return.list,
        objective = objective), return.list)
    } else {
      return.list = append(return.list, mean(marginals.child.1))
    }
    
    if (sd(marginals.child.2) >= 0.5 * abs(mean(marginals.child.1)) && length(marginals.child.2) > 2) {
      return.list = append(recursive_binary_split_marginal_effects(
        model = model, feature = feature, step.size = step.size,
        parent.node = child.node.2, return.list = return.list,
        objective = objective), return.list)
    } else {
      return.list = append(return.list, mean(marginals.child.2))
    }
  } else {
  }
  return(return.list)
}
