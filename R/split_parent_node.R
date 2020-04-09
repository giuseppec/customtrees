split_parent_node = function(Y, X, n.splits = 1, min.node.size = 10, objective, optimizer, ...) {
  assert_data_frame(X)
  #assert_choice(target.col, choices = colnames(data))
  assert_integerish(n.splits)
  assert_integerish(min.node.size)
  assert_function(objective, args = c("x", "y"))
  assert_function(optimizer, args = c("xval", "y"))

  # find best split points per feature
  opt.feature = lapply(X, function(feat) {
    optimizer(x = feat, y = Y, n.splits = n.splits, min.node.size = min.node.size, objective = objective, ...)
  })

  rbindlist(lapply(opt.feature, as.data.frame), idcol = "feature")
}
