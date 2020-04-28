# confidence_interval_sd_ratio = function(input.vector, target.ratio) {
# 
#   quantiles = quantile(input.vector, c(0.025, 0.975))
#   if (round(sd(input.vector), 10) == 0 || is.na(round(sd(input.vector), 10))) {
#     return(TRUE)
#   } else if (all((abs(quantiles) - mean(input.vector)) / sd(input.vector) < target.ratio)) {
#     return(TRUE)
#   } else {
#     return(FALSE)
#   }
# }

absolute_mean_sd_ratio = function(input.vector, target.ratio) {
  return(sd(input.vector) / abs(mean(input.vector)) <= target.ratio)
}