marginal_effects = function(model, data, feature, step.size) {
  prediction.before.intervention = predict(model, newdata = data)
  prediction.before.intervention = prediction.before.intervention$data$response
  data.modified = data
  data.modified[ , feature] = data.modified[ , feature] + step.size
  
  prediction.after.intervention = predict(model, newdata = data.modified)
  prediction.after.intervention = prediction.after.intervention$data$response
  
  marginal.effects = prediction.after.intervention - prediction.before.intervention
  return(marginal.effects)
}

