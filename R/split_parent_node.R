split_parent_node = function(Y, X, n.splits = 1, min.node.size = 1, optimizer,
  objective, ...) {
  assert_data_frame(X)
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

create_child_nodes = function(parent_node, parent_split_result) {
  
  split.feature = parent_split_result$feature[parent_split_result$best.split]
  if (length(split.feature) > 1) {
    # multiple features could provide the best split, select one randomly
    split.feature = sample(split.feature, 1)
  } else {
  }
  split.value = parent_split_result$split.points[which(parent_split_result$feature == split.feature)]
  
  child.1 = parent_node[parent_node[ , split.feature] <= split.value, ]
  child.2 = parent_node[parent_node[ , split.feature] > split.value, ]
  
  return(list(child.1, child.2))
}
