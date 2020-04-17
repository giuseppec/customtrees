marginal_effect_tree = function(model, feature, data, step.size, objective) {
  
  assign("return.list", list(), envir = .GlobalEnv)
  start.node = Node(id = "0", data = data)
  iter.count = 0
  assign("iter.count", iter.count, envir = .GlobalEnv)

  recursive_binary_split_marginal_effects(
    model = model, feature = feature, step.size = step.size,
    this.node = start.node, objective = objective)
}
# model = mod
# feature = "rm"
# step.size = 1
# objective = objective.function

recursive_binary_split_marginal_effects = function(model, feature, step.size,
                                                   this.node, objective) {
 
  marginals = marginal_effects(model, this.node$data, feature, step.size)
  this.node = set_node_ame(this.node, mean(marginals))
  this.node = set_node_n.observations(this.node, nrow(this.node$data))
  this.node = set_node_ame.sd(this.node, sd(marginals))
  
  iter.count = get("iter.count", envir = .GlobalEnv)
  iter.count = iter.count + 1
  print(paste("Recursive function call:", iter.count))
  assign("iter.count", iter.count, envir = .GlobalEnv)
  
  if (confidence_interval_sd_ratio(input.vector = marginals, target.ratio = 1)
      || nrow(this.node$data) <= 4) {
    # do not split node if 95% of the values within the node are located closer
    # to mean value than 1 standard deviation
    return.list = get("return.list", envir = .GlobalEnv)
    assign("return.list", append(return.list, list(this.node)), envir = .GlobalEnv)
    return(NULL)
  } else {
    # split node if within-node standard deviation too high
    split.result = split_parent_node(
      Y = marginals, X = this.node$data, objective = objective,
      optimizer = find_best_binary_split)
  
    if (nrow(split.result) == 0) {
      # node cannot be split
      return(NULL)
    }
    
    child.nodes = create_child_nodes(this.node, split.result)
    child.node.left = child.nodes[["child.1"]]
    child.node.right = child.nodes[["child.2"]]
    
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
        this.node = child.node.left, objective = objective)
    
    recursive_binary_split_marginal_effects(
      model = model, feature = feature, step.size = step.size,
      this.node = child.node.right, objective = objective)
    
  }
}
      
    # marginals.child.1 = marginal_effects(
  #   model, child.node.1$data, feature, step.size)
  # marginals.child.2 =  marginal_effects(
  #   model, child.node.2$data, feature, step.size)

  # return(list(this.node, recursive_binary_split_marginal_effects(
  #   # model = model, feature = feature, step.size = step.size,
  #   this.node = child.node.1, return.list = return.list,
  #   objective = objective)))
  # 
  

#   if (sd(marginals.child.1) >= 0.5 * abs(mean(marginals.child.1)) &&
#       length(marginals.child.1) > 2) {
#     return.list = append(return.list, recursive_binary_split_marginal_effects(
#       model = model, feature = feature, step.size = step.size,
#       this.node = child.node.1, return.list = return.list,
#       objective = objective))
#   } else {
#     
#     child.node.1 = set_node_ame(child.node.1, mean(marginals.child.1))
#     child.node.1 = set_node_n.observations(child.node.1, nrow(child.node.1$data))
#     return.list = append(return.list, list(child.node.1))
#   }
#   
#   if (sd(marginals.child.2) >= 0.5 * abs(mean(marginals.child.1)) &&
#       length(marginals.child.2) > 2) {
#     return.list = append(return.list, recursive_binary_split_marginal_effects(
#       model = model, feature = feature, step.size = step.size,
#       this.node = child.node.2, return.list = return.list,
#       objective = objective))
#   } else {
#     child.node.2 = set_node_ame(child.node.1, mean(marginals.child.2))
#     child.node.2 = set_node_n.observations(child.node.1, nrow(child.node.2$data))
#     return.list = append(return.list, list(child.node.2))
#   }
#   return(this.node)
# }
