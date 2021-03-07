# https://github.com/marjoleinF/gamtree#gam-based-recursive-partition-without-global-effects
library(ranger)
library(partykit)
library(mgcv)
data("BostonHousing", package = "mlbench")
BostonHousing$rad = as.factor(BostonHousing$rad)

mod = ranger(medv ~ ., data = BostonHousing, importance = "permutation")
mod_pred = predict(mod, data = BostonHousing)$predictions
sort(mod$variable.importance, decreasing = TRUE)


glmboost_node = function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  glmboost(y ~ x, ...)
}

glm_node = function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  glm(y ~ x + 0, start = start, ...)
}

step_node = function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  m = glm(y ~ 1, data = as.data.frame(x), start = start, ...)
  fac = grepl(paste(names(attr(x, "xlevels")), collapse = "|"), colnames(x))
  num_cols = colnames(x)[!fac]
  fac_cols = colnames(x)[fac]
  features = paste(c(num_cols, fac_cols), collapse = "+")
  form1 = as.formula(paste0("~", features,"+0"))
  
  m2 = step(m, direction = "forward", scope = list(lower=m, upper=form1))
}

gam_node = function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  fac = grepl(paste(names(attr(x, "xlevels")), collapse = "|"), colnames(x))
  num_cols = colnames(x)[!fac]
  fac_cols = colnames(x)[fac]
  features = paste(c(paste0("s(",num_cols, ")"), fac_cols), collapse = "+")
  form1 = as.formula(paste0("y~", features, "+0"))
  
  gam(form1, data = as.data.frame(x), ...)
  #gam(y ~ 0 + s(x), start = start, ...)
  #glm(y ~ x + poly(x, degree = 2, raw = TRUE), start = start, ...)
}

num_cols = colnames(sur_data)[sapply(sur_data, is.numeric)]
fac_cols = colnames(sur_data)[!sapply(sur_data, is.numeric)]
features = paste(c(paste0("poly(",num_cols, ", degree = 2, raw = TRUE)"), fac_cols), collapse = "+")
form_poly = as.formula(paste0("mod_pred~", features, "-1 | ", paste(c(num_cols, fac_cols), collapse = "+")))

features = paste(c(paste0("s(",num_cols, ")"), fac_cols), collapse = "+")
form_spline = as.formula(paste0("mod_pred~", features, "-1 | ", paste(c(num_cols, fac_cols), collapse = "+")))
form_interact = mod_pred ~ . + .:. - 1 | .

## set up a logistic regression tree
sur_data = subset(BostonHousing, select = -medv)
sur_mob = mob(mod_pred ~ . - 1 | ., data = sur_data, fit = glm_node)

sur_matrix = model.matrix(mod_pred ~ . - 1, data = sur_data, contrasts.arg = )
sur_mob_pred = predict(sur_mob, newdata = sur_data, type = "response")

mean((sur_mob_pred - mod_pred)^2)
plot(mod_pred, sur_mob_pred)

## visualization
partykit:::plot.modelparty(sur_mob)
partykit:::plot.modelparty(sur_mob, terminal_panel = node_bivplot)


## gamtree
# # install_github("marjoleinF/gamtree")
# num_cols = colnames(sur_data)[sapply(sur_data, is.numeric)]
# fac_cols = colnames(sur_data)[!sapply(sur_data, is.numeric)]
# features = paste(c(paste0("s(",num_cols, ")"), fac_cols), collapse = "+")
# form = as.formula(paste0("mod_pred~", features, "-1 | ", paste(c(num_cols, fac_cols), collapse = "+")))
# sur_gamtree = gamtree(form, data = cbind(mod_pred, sur_data), verbose = FALSE)
# sur_gamtree_pred = predict(sur_gamtree, sur_data)


## lm
sur_lm = lm(mod_pred ~ ., data = sur_data)
sur_lm_pred = predict(sur_lm, sur_data)

mean((sur_lm_pred - mod_pred)^2)
plot(mod_pred, sur_lm_pred)

## poly-lm
num_cols = colnames(sur_data)[sapply(sur_data, is.numeric)]
fac_cols = colnames(sur_data)[!sapply(sur_data, is.numeric)]
features = paste(c(paste0("poly(",num_cols, ", degree = 2, raw = FALSE)"), fac_cols), collapse = "+")
form = as.formula(paste0("mod_pred~", features))
sur_plm = lm(form, data = sur_data)
sur_plm_pred = predict(sur_plm, sur_data)

mean((sur_plm_pred - mod_pred)^2)
plot(mod_pred, sur_plm_pred)

## gam
library(mgcv)
num_cols = colnames(sur_data)[sapply(sur_data, is.numeric)]
fac_cols = colnames(sur_data)[!sapply(sur_data, is.numeric)]
features = paste(c(paste0("s(",num_cols, ")"), fac_cols), collapse = "+")
form = as.formula(paste0("mod_pred~", features))

sur_gam = gam(form, data = sur_data)
sur_gam_pred = predict(sur_gam, sur_data)

mean((sur_gam_pred - mod_pred)^2)
plot(mod_pred, sur_gam_pred)



data("PimaIndiansDiabetes", package = "mlbench")

## a simple basic fitting function (of type 1) for a logistic regression
logit <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  glm(y ~ 0 + x, family = binomial, start = start, ...)
}

## set up a logistic regression tree
pid_tree <- mob(diabetes ~ glucose | pregnant + pressure + triceps + insulin +
    mass + pedigree + age, data = PimaIndiansDiabetes, fit = logit)
## see lmtree() and glmtree() for interfaces with more efficient fitting functions

## print tree
plot(pid_tree)

plot.modelparty(pid_tree, terminal_panel = node_bivplot)
plot.modelparty(pid_tree, terminal_panel = node_scatterplot)

FUN <- function(x) {
  cf <- x$coefficients
  cf <- matrix(cf, ncol = 1, dimnames = list(names(cf), 
    ""))
  c(sprintf("n = %s", x$nobs), "Estimated parameters:", 
    strwrap(capture.output(print(cf, digits = 4L))[-1L]))
}
node = do.call("node_terminal", c(list(obj = pid_tree, FUN = FUN)))
plot.party(pid_tree, terminal_panel = node)






























# https://github.com/marjoleinF/gamtree#gam-based-recursive-partition-without-global-effects
library(ranger)
library(partykit)
library(mgcv)
data("BostonHousing", package = "mlbench")
BostonHousing$rad = as.factor(BostonHousing$rad)

mod = ranger(medv ~ ., data = BostonHousing, importance = "permutation")
mod_pred = predict(mod, data = BostonHousing)$predictions
sort(mod$variable.importance, decreasing = TRUE)

step_node = function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  m = glm(y ~ 1, start = start, ...)
  m2 = step(m, direction = "forward", scope = list(lower=m, upper=~x+0))
}

glm_node = function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  glm(y ~ x + 0, start = start, ...)
}

model_node = function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  # #gam(y ~ 0 + s(x), start = start, ...)
  # #glm(y ~ x + poly(x, degree = 2, raw = FALSE), start = start, ...)
  # fac = grepl(paste(names(attr(x, "xlevels")), collapse = "|"), colnames(x))
  # num_cols = colnames(x)[!fac]
  # fac_cols = colnames(x)[fac]
  # features = paste(c(paste0("s(", num_cols, ")"), fac_cols), collapse = "+")
  # form = as.formula(paste0("y~", features))
  # 
  # gam(form, data = as.data.frame(x), ...)
  gam(y~x, ...)
}

num_cols = colnames(sur_data)[sapply(sur_data, is.numeric)]
fac_cols = colnames(sur_data)[!sapply(sur_data, is.numeric)]
features = paste(c(paste0("poly(",num_cols, ", degree = 2, raw = FALSE)"), fac_cols), collapse = "+")
form = as.formula(paste0("mod_pred~", features, "-1 | ", paste(c(num_cols, fac_cols), collapse = "+")))
form = mod_pred ~ . + .:. - 1 | .

## set up a logistic regression tree
sur_data = subset(BostonHousing, select = -medv)
sur_mob = mob(form, data = sur_data, fit = glm_node)#, control = mob_control(bonferroni = FALSE))
sur_mob_pred = predict(sur_mob, data = sur_data, type = "response")

mean((sur_mob_pred - mod_pred)^2)
plot(mod_pred, sur_mob_pred)

## visualization
partykit:::plot.modelparty(sur_mob)
partykit:::plot.modelparty(sur_mob, terminal_panel = node_bivplot)


## gamtree
# # install_github("marjoleinF/gamtree")
# num_cols = colnames(sur_data)[sapply(sur_data, is.numeric)]
# fac_cols = colnames(sur_data)[!sapply(sur_data, is.numeric)]
# features = paste(c(paste0("s(",num_cols, ")"), fac_cols), collapse = "+")
# form = as.formula(paste0("mod_pred~", features, "-1 | ", paste(c(num_cols, fac_cols), collapse = "+")))
# sur_gamtree = gamtree(form, data = cbind(mod_pred, sur_data), verbose = FALSE)
# sur_gamtree_pred = predict(sur_gamtree, sur_data)


## lm
sur_lm = lm(mod_pred ~ ., data = sur_data)
sur_lm_pred = predict(sur_lm, sur_data)

mean((sur_lm_pred - mod_pred)^2)
plot(mod_pred, sur_lm_pred)

## poly-lm
num_cols = colnames(sur_data)[sapply(sur_data, is.numeric)]
fac_cols = colnames(sur_data)[!sapply(sur_data, is.numeric)]
features = paste(c(paste0("poly(",num_cols, ", degree = 2, raw = FALSE)"), fac_cols), collapse = "+")
form = as.formula(paste0("mod_pred~", features))
sur_plm = lm(form, data = sur_data)
sur_plm_pred = predict(sur_plm, sur_data)

mean((sur_plm_pred - mod_pred)^2)
plot(mod_pred, sur_plm_pred)

## gam
library(mgcv)
num_cols = colnames(sur_data)[sapply(sur_data, is.numeric)]
fac_cols = colnames(sur_data)[!sapply(sur_data, is.numeric)]
features = paste(c(paste0("s(",num_cols, ")"), fac_cols), collapse = "+")
form = as.formula(paste0("mod_pred~", features))

sur_gam = gam(form, data = sur_data)
sur_gam_pred = predict(sur_gam, sur_data)

mean((sur_gam_pred - mod_pred)^2)
plot(mod_pred, sur_gam_pred)



data("PimaIndiansDiabetes", package = "mlbench")

## a simple basic fitting function (of type 1) for a logistic regression
logit <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  glm(y ~ 0 + x, family = binomial, start = start, ...)
}

## set up a logistic regression tree
pid_tree <- mob(diabetes ~ glucose | pregnant + pressure + triceps + insulin +
    mass + pedigree + age, data = PimaIndiansDiabetes, fit = logit)
## see lmtree() and glmtree() for interfaces with more efficient fitting functions

## print tree
plot(pid_tree)

plot.modelparty(pid_tree, terminal_panel = node_bivplot)
plot.modelparty(pid_tree, terminal_panel = node_scatterplot)

FUN <- function(x) {
  cf <- x$coefficients
  cf <- matrix(cf, ncol = 1, dimnames = list(names(cf), 
    ""))
  c(sprintf("n = %s", x$nobs), "Estimated parameters:", 
    strwrap(capture.output(print(cf, digits = 4L))[-1L]))
}
node = do.call("node_terminal", c(list(obj = pid_tree, FUN = FUN)))
plot.party(pid_tree, terminal_panel = node)
