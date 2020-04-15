Node = function(id = NULL, data, split.feature = NULL, split.value = NULL,
                child.left = NULL, child.right = NULL, ame = NULL) {
  node.values = list(
    "id" = id, "data" = data, "split.feature" = split.feature,
    "split.value" = split.value, "child.left" = child.left,
    "child.right" = child.right, "ame" = ame)
  attr(node.values, "class") = "Node"
  return(node.values)
}

get_node_id = function(Node.obj) {
  Node.obj$id
}

get_node_data = function(Node.obj) {
  Node.obj$data
}

get_node_child.left = function(Node.obj) {
  Node.obj$child.left
}

get_node_child.right = function(Node.obj) {
  Node.obj$child.right
}

get_node_split.feature = function(Node.obj) {
  Node.obj$split.feature
}

get_node_split.value = function(Node.obj) {
  Node.obj$split.value
}

get_node_ame = function(Node.obj) {
  Node.obj$ame
}

set_node_id = function(Node.obj, note.id) {
  Node.obj$id = note.id
  return(Node.obj)
}

set_node_data = function(Node.obj, note.data) {
  Node.obj$data = note.data
  return(Node.obj)
}
set_node_child.left = function(Node.obj, child.left) {
  Node.obj$child.left = child.left
  return(Node.obj)
}

set_node_child.right = function(Node.obj, child.right) {
  Node.obj$child.right = child.right
  return(Node.obj)
}

set_node_split.feature = function(Node.obj, split.feature) {
  Node.obj$split.feature = split.feature
  return(Node.obj)
}

set_node_split.value = function(Node.obj, split.value) {
  Node.obj$split.value = split.value
  return(Node.obj)
}

set_node_ame = function(Node.obj, ame) {
  Node.obj$ame = ame
  return(Node.obj)
}
