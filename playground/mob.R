# https://github.com/marjoleinF/gamtree#gam-based-recursive-partition-without-global-effects
library(ranger)
library(partykit)
library(mgcv)
library(leaps)
library(glmnet)
data("BostonHousing", package = "mlbench")
BostonHousing$rad = as.factor(BostonHousing$rad)

# fit black-box model
mod = ranger(medv ~ ., data = BostonHousing, importance = "permutation")
mod_pred = predict(mod, data = BostonHousing)$predictions
sort(mod$variable.importance, decreasing = TRUE)

glm_node = function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  glm(y ~ x + 0, start = start, ...)
}

leaps_node = function(y, x, start = NULL, weights = NULL, offset = NULL, ...){
  models = regsubsets(x = x, y = y, nvmax = 4, really.big = TRUE)
  m = summary(models)
  id = which.max(m$adjr2)
  feats = m$which[id,-1]
  x = x[,feats]
  glm(y ~ x + 0, start = start, ...)
}

lasso_node = function(y, x, start = NULL, weights = NULL, offset = NULL, ...){
  mod = cv.glmnet(x = x, y = y, alpha = 1, intercept = FALSE)
  coef = coef(mod, s = "lambda.1se")
  feats = row.names(coef)[as.vector(coef!=0)]
  x = x[,feats]
  glm(y ~ x + 0, start = start, ...)
}

# step_node = function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
#   m = glm(y ~ 1, data = as.data.frame(x), start = start, ...)
#   fac = grepl(paste(names(attr(x, "xlevels")), collapse = "|"), colnames(x))
#   num_cols = colnames(x)[!fac]
#   fac_cols = colnames(x)[fac]
#   features = paste(c(num_cols, fac_cols), collapse = "+")
#   form1 = as.formula(paste0("~", features,"+0"))
# 
#   m2 = step(m, direction = "forward", scope = list(lower=m, upper=form1), steps = 3)
# }

# non-linear effects for numeric features
sur_data = subset(BostonHousing, select = -medv)
num_cols = colnames(sur_data)[sapply(sur_data, is.numeric)]
fac_cols = colnames(sur_data)[!sapply(sur_data, is.numeric)]
feat_poly = paste(c(paste0("poly(",num_cols, ", degree = 2, raw = TRUE)"), fac_cols), collapse = "+")
feat_spline = paste(c(paste0("s(",num_cols, ")"), fac_cols), collapse = "+")
feat_sq = paste(c(num_cols, paste0("I(",num_cols, "^2)"), fac_cols), collapse = "+")

# crate some different formulas for the use in mob()
form_poly = as.formula(paste0("mod_pred~", feat_poly, "-1 | ", paste(c(num_cols, fac_cols), collapse = "+")))
form_spline = as.formula(paste0("mod_pred~", feat_spline, "-1 | ", paste(c(num_cols, fac_cols), collapse = "+")))
form_sq = as.formula(paste0("mod_pred~", feat_sq, "-1 | ", paste(c(num_cols, fac_cols), collapse = "+")))
form_interact = mod_pred ~ . + .:. - 1 | .
form = mod_pred ~ . - 1 | .

## 0. fit mob surrogate model
##############################
# max of 10% of the data or 4 times number of parameters
npar_formula = as.formula(paste0("~",gsub("[|].*", "", as.character(form)[3])))
npars = ncol(model.matrix(npar_formula, sur_data))
minsize = 3*npars #max(c(5*npars, floor(0.1*nrow(sur_data))))

sur_mob = mob(form, data = cbind(mod_pred, sur_data), fit = glm_node, 
  control = mob_control(alpha = 0.2, minsize = minsize)) 
sur_mob_pred = predict(sur_mob, newdata = sur_data, type = "response")
# node sizes:
table(predict(sur_mob, type = "node"))

# check approx. of mob surrogate model
mean((sur_mob_pred - mod_pred)^2)
plot(mod_pred, sur_mob_pred)

# visualization of mob surrogate model
partykit:::plot.modelparty(sur_mob)
partykit:::plot.modelparty(sur_mob, terminal_panel = node_bivplot) # buggy

nodes = nodeapply(sur_mob, ids = nodeids(sur_mob, terminal = TRUE), info_node)
names(nodes)
# 
# library(effects)
# eff.pres = allEffects(nodes[[2]]$object)
# 
# plot(eff.pres)
# nodes[[2]]

## 1. fit CART surrogate model
#############################
library(rpart)
sur_rp = rpart(mod_pred ~ ., data = sur_data)
sur_rp_pred = predict(sur_rp, sur_data)

mean((sur_rp_pred - mod_pred)^2)
plot(mod_pred, sur_rp_pred)

## 2. fit lm surrogate model
#############################
sur_lm = lm(mod_pred ~ ., data = sur_data)
sur_lm_pred = predict(sur_lm, sur_data)

mean((sur_lm_pred - mod_pred)^2)
plot(mod_pred, sur_lm_pred)

## 3. fit poly-lm surrogate model
##################################
num_cols = colnames(sur_data)[sapply(sur_data, is.numeric)]
fac_cols = colnames(sur_data)[!sapply(sur_data, is.numeric)]
features = paste(c(paste0("poly(",num_cols, ", degree = 2, raw = FALSE)"), fac_cols), collapse = "+")
form = as.formula(paste0("mod_pred~", features))
sur_plm = lm(form, data = sur_data)
sur_plm_pred = predict(sur_plm, sur_data)

mean((sur_plm_pred - mod_pred)^2)
plot(mod_pred, sur_plm_pred)

## 4. fit gam surrogate model
##################################
library(mgcv)
num_cols = colnames(sur_data)[sapply(sur_data, is.numeric)]
fac_cols = colnames(sur_data)[!sapply(sur_data, is.numeric)]
features = paste(c(paste0("s(",num_cols, ")"), fac_cols), collapse = "+")
form = as.formula(paste0("mod_pred~", features))
sur_gam = gam(form, data = sur_data)
sur_gam_pred = predict(sur_gam, sur_data)

mean((sur_gam_pred - mod_pred)^2)
plot(mod_pred, sur_gam_pred)
