Node = function(id = NULL, depth = NULL, data, split.feature = NULL, split.value = NULL,
                parent.node = NULL, child.left = NULL, child.right = NULL, ame = NULL,
                ame.sd = NULL, n.observations = NULL) {
  node.values = list(
    "id" = id, "depth" = depth, "data" = data, "split.feature" = split.feature,
    "split.value" = split.value, "parent.node" = parent.node, "child.left" = child.left,
    "child.right" = child.right, "ame" = ame, "ame.sd" = ame.sd,
    "n.observations" = n.observations)
  attr(node.values, "class") = "Node"
  return(node.values)
}

get_node_id = function(Node.obj) {
  Node.obj$id
}

get_node_depth = function(Node.obj) {
  Node.obj$depth
}

get_node_data = function(Node.obj) {
  Node.obj$data
}

get_node_parent.node = function(Node.obj) {
  Node.obj$parent.node
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

get_node_n.observations = function(Node.obj) {
  Node.obj$n.observations
}

get_node_ame.sd = function(Node.obj) {
  Node.obj$ame.sd
}

set_node_id = function(Node.obj, node.id) {
  Node.obj$id = node.id
  return(Node.obj)
}

set_node_depth = function(Node.obj, node.depth) {
  Node.obj$depth = node.depth
  return(Node.obj)
}

set_node_data = function(Node.obj, node.data) {
  Node.obj$data = node.data
  return(Node.obj)
}

set_node_parent.node = function(Node.obj, parent.node) {
  Node.obj$parent.node = parent.node
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

set_node_n.observations = function(Node.obj, n.observations) {
  Node.obj$n.observations = n.observations
  return(Node.obj)
}

set_node_ame.sd = function(Node.obj, ame.sd) {
  Node.obj$ame.sd = ame.sd
  return(Node.obj)
}
