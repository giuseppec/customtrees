split_parent_node = function(Y, X, n.splits = 1, min.node.size = 1, optimizer, dist.between = "fre", lambda = NULL,
  objective, ...) {
  assert_data_frame(X)
  #assert_choice(target.col, choices = colnames(data))
  assert_integerish(n.splits)
  assert_integerish(min.node.size)
  assert_function(objective, args = c("y", "x", "requires.x"))
  assert_function(optimizer, args = c("xval", "y"))

  # find best split points per feature
  opt.feature = lapply(X, function(feat) {
    optimizer(x = feat, y = Y, n.splits = n.splits, min.node.size = min.node.size, dist.between = dist.between, lambda = lambda,
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
