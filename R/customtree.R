customtree = function(data, target.cols, n.splits = 1, min.node.size = 10, objective, optimizer) {
  assert_data_frame(data)
  assert_choice(target.col, choices = colnames(data))
  assert_integerish(n.splits)
  assert_integerish(min.node.size)
  assert_function(objective, args = c("x", "y"))
  assert_function(optimizer)

}
