marginal_effects = function(model, data, feature, step.size) {
  prediction.before.intervention = predict(model$learner.model, data)
  data.modified = data
  data.modified[ , feature] = data.modified[ , feature] + step.size
  
  prediction.after.intervention = predict(model$learner.model, data.modified)
  
  marginal.effects = prediction.after.intervention - prediction.before.intervention
  return(marginal.effects)
}