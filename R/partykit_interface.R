
party_interface = function(tree.summary, data) {
  
  party.node.list = list()
  depth = max(tree.summary$depth)
  nodes = tree.summary[which(tree.summary$depth == depth), ]
  nodes.id = as.numeric(as.character(nodes$id))
  party.nodes = lapply(nodes.id, FUN = function(node) {
    partynode(node)
  })
  party.node.list = append(party.node.list, party.nodes)
  depth = depth - 1
  while (depth > 0) {
    
    nodes = tree.summary[tree.summary$depth == depth, ]
    parent.nodes = nodes[!is.na(nodes$split.feature), ]
    leaf.nodes = nodes[is.na(nodes$split.feature), ]
    
    parent.party.nodes = lapply(1:nrow(parent.nodes), FUN = function(i) {
      node = parent.nodes[i, ]
      node.id = as.character(node$id)
      node.split.feature = node$split.feature
      node.split.value = node$split.value
      node.split.feature.index = which(features == node.split.feature)
      
      child.party.nodes = lapply(party.node.list, FUN = function(x) {
        
        x.splitstring = unlist(strsplit(as.character(x$id), ""))
        x.parent = paste0(x.splitstring[-length(x.splitstring)], collapse = "")
        if (node.id == x.parent) {
          return(x)
        } else {
        }})
      
      child.party.nodes = child.party.nodes[!unlist(lapply(child.party.nodes, is.null))]
      # drop NULL values
      
      party.node = partynode(
        as.numeric(node.id), split = partysplit(node.split.feature.index, node.split.value),
        kids = child.party.nodes)
      return(party.node)
    })
    
    if (nrow(leaf.nodes) > 0) {
      leaf.party.nodes = lapply(1:nrow(leaf.nodes), FUN = function(i) {
        node = leaf.nodes[i, ]
        node.id = as.numeric(as.character(node$id))
        party.node = partynode(
          node.id)
        return(party.node)
      })
    } else {
    }
    party.node.list = append(party.node.list, parent.party.nodes)
    party.node.list = append(party.node.list, leaf.party.nodes)
    
    depth = depth - 1
  }
  
  tree.index = which(unlist(lapply(party.node.list, FUN = function(x) {
    x$id == 1
  })))
  
  party.tree = party.node.list[[tree.index]]
  party.tree = party(party.tree, data)
  return(party.tree)
}