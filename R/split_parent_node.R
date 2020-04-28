split_parent_node = function(Y, X, n.splits = 1, min.node.size = 1, optimizer,
  objective, ...) {
  # assert_data_frame(X)
  #assert_choice(target.col, choices = colnames(data))
  assert_integerish(n.splits)
  assert_integerish(min.node.size)
  assert_function(objective, args = c("y", "x", "requires.x"))
  assert_function(optimizer, args = c("xval", "y"))

  # find best split points per feature
  opt.feature = lapply(X, function(feat) {
    optimizer(x = feat, y = Y, n.splits = n.splits, min.node.size = min.node.size,
      objective = objective, ...)
  })

  result = rbindlist(lapply(opt.feature, as.data.frame), idcol = "feature")
  result = result[, .(split.points = list(split.points)), by = c("feature", "objective.value"), with = TRUE]
  result$best.split = result$objective.value == min(result$objective.value)
  #result = result[, best.split := objective.value == min(objective.value)]
  return(result)
}

generate_node_index = function(Y, X, result) {
  assert_data_table(result)
  # TODO: fix bug if more than one feature have the same best objective
  feature = unique(result$feature[result$best.split])
  split.points = unlist(result$split.points[result$best.split])

  if (is.vector(X))
    xval = X else
      xval = X[, feature]

  cuts = c(min(xval), split.points, max(xval))
  sp = cut(xval, breaks = unique(cuts), include.lowest = TRUE)
  #levels(sp) = paste0(feature, " in ", levels(sp))

  return(list(class = sp, index = split(seq_along(xval), sp)))
}

create_child_nodes = function(parent.node, parent.split.result) {
  
  split.data = get_split_data(parent.split.result)
  split.feature = split.data[["split.feature"]]
  split.value = split.data[["split.value"]]
  
  parent.node.data = get_node_data(parent.node)
  parent.node.id = get_node_id(parent.node)
  
  child.1.data = parent.node.data[parent.node.data[ , split.feature] <= split.value, ]
  child.1.id = paste0(as.character(parent.node.id), 0)
  child.1.depth = unlist(lapply(strsplit(child.1.id, ""), FUN = length))
  
  child.2.data = parent.node.data[parent.node.data[ , split.feature] > split.value, ]
  child.2.id = paste0(as.character(parent.node.id), 1)
  child.2.depth = unlist(lapply(strsplit(child.2.id, ""), FUN = length))
  
  child.1 = Node(
    id = child.1.id,
    depth = child.1.depth,
    data = child.1.data)
  
  child.2 = Node(
    id = child.2.id,
    depth = child.2.depth,
    data = child.2.data)
  
  return(list("child.1" = child.1, "child.2" = child.2))
}

get_split_data = function(parent.split.result) {
  
  split.feature = parent.split.result$feature[parent.split.result$best.split]
  if (length(split.feature) > 1) {
    # multiple features could provide the best split, select one randomly
    split.feature = sample(split.feature, 1)
  } else {
    # if only one feature provides the best split, do nothing
  }
  split.value = parent.split.result$split.points[[
    which(parent.split.result$feature == split.feature)]]
  objective.value = parent.split.result$objective.value[
    which(parent.split.result$feature == split.feature)]
  
  return(list("split.feature" = split.feature, "split.value" = split.value,
              "objective.value" = objective.value))
}
