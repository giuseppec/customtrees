# interaction effects for one feature
interactions_per_feature = function(x, y, feat, n.splits, min.node.size, optimizer, objective.split, objective.total, improve.first.split, improve.n.splits, ...){
  res = vector("list", length = n.splits)
  for (i in 1:n.splits) {
    res[[i]] = split_parent_node(Y = y, X = x[,-which(colnames(x) == feat)], n.splits = i, min.node.size = min.node.size, 
      optimizer = optimizer, objective = objective.split, feat = feat, ...)
    if (i == 1) {#is.na(unlist(improve)) in case improve gets NaN (0/0)
      obj.total = objective.total(y = y, xval = x, feat = feat, ...)
      improve = (obj.total - res[[i]][which(res[[i]]$best.split == TRUE), "objective.value"])/obj.total
      if (improve < improve.first.split | is.nan(unlist(improve))) return(NULL)
    }
    if (i > 1) {
      obj.sub1 = res[[i - 1]][which(res[[i - 1]]$best.split == TRUE), "objective.value"]
      obj.sub2 = res[[i]][which(res[[i - 1]]$best.split == TRUE), "objective.value"]
      improve = (obj.sub1 - obj.sub2)/(obj.sub1)
      if (improve < improve.n.splits | is.nan(unlist(improve))) return(list(res[(i - 1)], obj.total))
    }
  }
  return(list(res[[i]], obj.total))
}


# find potential interactions between all features in X 
find_potential_interactions = function(x, model, optimizer, objective.split, objective.total, n.splits = 3, min.node.size = 1, improve.first.split = 0.1, improve.n.splits = 0.1, ...){
  p = ncol(x)
  for (i in 1:p) {
    effect = FeatureEffect$new(model, method = "ice", grid.size = 20, feature = colnames(x)[i])
    # Get ICE values and arrange them in a horizontal matrix
    colnames(effect$results)[1] = ".borders"
    Y = spread(effect$results, .borders, .value)
    Y = Y[, setdiff(colnames(Y), c(".type", ".id"))]
    for (j in 1:nrow(Y)) {
      Y[j,] = as.numeric(unname(Y[j,])) - mean(as.numeric(unname(Y[j,])))
    }
    result.splits = interactions_per_feature(x = x, y = Y, feat = colnames(x)[i], n.splits = n.splits, min.node.size = min.node.size,
      optimizer = optimizer, objective.split = objective.split, objective.total = objective.total, 
      improve.first.split = improve.first.split, improve.n.splits = improve.n.splits, ...)
    
    #if(i == 1 & is.null(result.splits)) next
    if (!is.null(result.splits)) {
      result.sub = data.frame(result.splits[[1]])
      result.sub$improvement = (result.splits[[2]] - result.sub$objective.value)/result.splits[[2]]
      result.sub$ICEfeature = colnames(x)[i]
      if (exists("res") == FALSE) res = result.sub
      else res = bind_rows(res, result.sub)
    }
  }
  return(res)
}


# find "true" interactions 
find_true_interactions = function(x, model, optimizer, objective.split = SS_fre_filtered, objective.total = SS_fre, n.splits = 3, min.node.size = 10, improve.first.split = 0.1, improve.n.splits = 0.1, interaction.factor = 0.1, ...){
   # get potential interactions
  potential.interactions = find_potential_interactions(x = x, model = model, optimizer = optimizer, objective.split = objective.split, objective.total = objective.total, 
    n.splits = n.splits, min.node.size = min.node.size, ...)
  # user only interactions that have a higher effect than interaction.factor and that are two-sided
  true.interactions = potential.interactions[which(potential.interactions$improvement > interaction.factor),]
  
  ind = integer(length = 0)
  for (i in 1:nrow(true.interactions)) {
    n.row = nrow(true.interactions[which(which(true.interactions$feature %in% true.interactions$ICEfeature[i]) %in% which(true.interactions$ICEfeature %in% true.interactions$feature[i])),])
    if (n.row == 0) {
      ind = c(ind, i)
    }
  }
  if (length(ind) != 0) true.interactions = true.interactions[-ind,]
  
  return(true.interactions[order(true.interactions$improvement, decreasing = TRUE),])
}

