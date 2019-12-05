data(iris)

Y = as.numeric(iris$Species)
X = as.matrix(iris[, -c(5)])
X_stumps = make_stumps_matrix(X, T)


learner_iris = VEBBoostMultiClassLearner$new()
learner_iris$X = X_stumps
learner_iris$Y = Y
learner_iris$learners = list(
  learner1 = VEBBoostNode$new("mu_0", fitFunction = fitFn, predFunction = predFn, currentFit = list(mu1 = 0, mu2 = 1, KL_div = 1), family = 'binomial'), 
  learner2 = VEBBoostNode$new("mu_0", fitFunction = fitFn, predFunction = predFn, currentFit = list(mu1 = 0, mu2 = 1, KL_div = 1), family = 'binomial'), 
  learner3 = VEBBoostNode$new("mu_0", fitFunction = fitFn, predFunction = predFn, currentFit = list(mu1 = 0, mu2 = 1, KL_div = 1), family = 'binomial')
)

learner_iris$learners$learner1$ensemble = learner_iris
learner_iris$learners$learner2$ensemble = learner_iris
learner_iris$learners$learner3$ensemble = learner_iris

learner_iris$learners$learner1$Y = as.numeric(Y == 1)
learner_iris$learners$learner2$Y = as.numeric(Y == 2)
learner_iris$learners$learner3$Y = as.numeric(Y == 3)

learner_iris$convergeFit()

learner_iris$convergeFitAll(tol = 1e-1)


while ((abs(tail(tail(learner_iris$ELBO_progress, 1)[[1]], 1) - tail(tail(learner_iris$ELBO_progress, 2)[[1]], 1)) > tol) && (sum(sapply(learner_iris$learners, function(x) length(Traverse(x, filterFun = function(node) node$isLeaf & !node$isLocked)))) > 0)) {
  learner_iris$convergeFitAll(1e-1, V_tol = 1e-3)
}
