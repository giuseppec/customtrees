record_split_sequence = function(tree, leaf.id) {
  
  id.split = unlist(strsplit(x = leaf.id, split = "", fixed = TRUE))
  split.list = list()
  id = character(1)
  
  for (i in 1:(length(id.split) - 1)) {
    
    id = paste0(id, id.split[i])
    split = tree[tree$id == id, c("split.feature", "split.value")]
    if (id.split[i + 1] == "0") {
      split$direction = "left" 
    } else {
      split$direction = "right"  
    }
    split.list = append(split.list, list(split))
  }
  split.sequence = data.frame(do.call(rbind, split.list))
  split.sequence$split.feature = as.character(split.sequence$split.feature)
  split.sequence$split.value = as.numeric(as.character(split.sequence$split.value))
  return(split.sequence)
}



get_feature_subspace = function(data, split.sequence) {
  
  split.features = unique(split.sequence$split.feature)
  
  feature.subspace.list = lapply(split.features, FUN = function(feature) {
    feature.splits = split.sequence[split.sequence$split.feature == feature, ]
    
    if ("left" %in% feature.splits$direction) {
      boundary.right = min(feature.splits$split.value[feature.splits$direction == "left"])
    } else {
      boundary.right = max(data[[feature]])
    }
    if ("right" %in% feature.splits$direction) {
      boundary.left = max(feature.splits$split.value[feature.splits$direction == "right"])
    } else {
      boundary.left = min(data[[feature]])
    }
    feature.subspace = paste0("[", boundary.left, ",", boundary.right, "]")
    return(feature.subspace)
  })
  names(feature.subspace.list) = split.features
  feature.subspace.data.frame = do.call(cbind.data.frame, feature.subspace.list)
  return(feature.subspace.data.frame)
}

get_global_subspace = function(data, leaf.nodes) {
  
  feature.subspaces = lapply(leaf.nodes$id, FUN = function(id) {
    split.sequence = record_split_sequence(tree = tree.summary, leaf.id = id)
    feature.subspace.list = get_feature_subspace(data, split.sequence)
    return(feature.subspace.list)
  })
  
  feature.subspaces = lapply(feature.subspaces, FUN = function(feature.subspace) {
    complementary.features = colnames(data)[!colnames(data) %in% colnames(feature.subspace)]
    feature.subspace[complementary.features] = NA 
    return(feature.subspace)
  })
  global.subspace = do.call(rbind.data.frame, feature.subspaces)
  global.subspace = apply(global.subspace, MARGIN = 2, FUN = as.character)
  
  for (j in 1:ncol(global.subspace)) {
    
    feature = colnames(global.subspace)[j]
    feature.range = range(data[[feature]])
    feature.range =  paste0("[", feature.range[1], ",", feature.range[2], "]")
    for (i in 1:nrow(global.subspace)) {
      obs = global.subspace[i, feature]
      if (is.na(obs)) {
        global.subspace[i, feature] = feature.range
      } else {
      }
    }
  }
  
  global.subspace = as.data.frame(global.subspace)
  global.subspace = cbind.data.frame(leaf.nodes, global.subspace)
  
  return(global.subspace)
}
