marginal_effect_tree = function(model, feature, data, step.size, objective, target.ratio.sd.mean,
                                min.node.size) {
  
  assign("return.list", list(), envir = .GlobalEnv)
  start.node = Node(id = "1", "parent.node" = NA, "depth" = 1, data = data)
  iter.count = 0
  assign("iter.count", iter.count, envir = .GlobalEnv)

  recursive_binary_split_marginal_effects(
    model = model, feature = feature, step.size = step.size,
    this.node = start.node, objective = objective, target.ratio.sd.mean = target.ratio.sd.mean,
    min.node.size = min.node.size)
}

recursive_binary_split_marginal_effects = function(model, feature, step.size,
                                                   this.node, objective, target.ratio.sd.mean,
                                                   min.node.size) {
 
  marginals = marginal_effects(model, this.node$data, feature, step.size)
  this.node = set_node_ame(this.node, mean(marginals))
  this.node = set_node_n.observations(this.node, nrow(this.node$data))
  this.node = set_node_ame.sd(this.node, sd(marginals))
  this.node = set_node_split.feature(this.node, NA)
  this.node = set_node_split.value(this.node, NA)
  this.node = set_node_child.left(this.node, NA)
  this.node = set_node_child.right(this.node, NA)
  
  iter.count = get("iter.count", envir = .GlobalEnv)
  iter.count = iter.count + 1
  print(paste("Recursive function call:", iter.count))
  assign("iter.count", iter.count, envir = .GlobalEnv)
  
  if (absolute_mean_sd_ratio(input.vector = marginals, target.ratio = target.ratio.sd.mean)) {
    return.list = get("return.list", envir = .GlobalEnv)
    assign(
      "return.list", append(return.list, list(this.node)), envir = .GlobalEnv)
    return(NULL)
  } else {
    print(nrow(this.node$data))
    # split node if within-node standard deviation too high
    split.result = split_parent_node(
      Y = marginals, X = this.node$data, objective = objective,
      optimizer = find_best_binary_split, min.node.size = min.node.size)
  }
  
  if (nrow(split.result) == 0) {
    # node cannot be split
    return(NULL)
  } else if (get_split_data(split.result)[["objective.value"]] == Inf) {
    return.list = get("return.list", envir = .GlobalEnv)
    assign(
      "return.list", append(return.list, list(this.node)), envir = .GlobalEnv)
    return(NULL)
  }
  
  child.nodes = create_child_nodes(this.node, split.result)
  child.node.left = child.nodes[["child.left"]]
  child.node.left = set_node_parent.node(child.node.left, parent.node = this.node$id)
  child.node.right = child.nodes[["child.right"]]
  child.node.right = set_node_parent.node(child.node.right, parent.node = this.node$id)
  
  this.node = set_node_split.feature(
    this.node, get_split_data(split.result)[["split.feature"]])
  this.node = set_node_split.value(
    this.node, get_split_data(split.result)[["split.value"]])
  this.node = set_node_child.left(this.node, child.node.left$id)
  this.node = set_node_child.right(this.node, child.node.right$id)
  
  return.list = get("return.list", envir = .GlobalEnv)
  assign("return.list", append(return.list, list(this.node)), envir = .GlobalEnv)
  
  recursive_binary_split_marginal_effects(
    model = model, feature = feature, step.size = step.size,
    this.node = child.node.left, objective = objective,
    target.ratio.sd.mean = target.ratio.sd.mean, min.node.size = min.node.size)
  
  recursive_binary_split_marginal_effects(
    model = model, feature = feature, step.size = step.size,
    this.node = child.node.right, objective = objective,
    target.ratio.sd.mean = target.ratio.sd.mean, min.node.size = min.node.size)
}

