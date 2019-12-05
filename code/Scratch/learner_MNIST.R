train <- read.csv("D:/School/Fall 2019/Matthew Stephens Research/Kaggle MNIST Data/train.csv", stringsAsFactors=FALSE)


Y = train$label
X = as.matrix(train[, -1])
X = X / 255
X_stumps = make_stumps_matrix(X, F)


learner_MNIST = VEBBoostMultiClassLearner$new()
learner_MNIST$X = X_stumps
learner_MNIST$Y = Y
learner_MNIST$learners = list(
  learner0 = VEBBoostNode$new("mu_0", fitFunction = fitFn, predFunction = predFn, currentFit = list(mu1 = 0, mu2 = 1, KL_div = 1), family = 'binomial'), 
  learner1 = VEBBoostNode$new("mu_0", fitFunction = fitFn, predFunction = predFn, currentFit = list(mu1 = 0, mu2 = 1, KL_div = 1), family = 'binomial'), 
  learner2 = VEBBoostNode$new("mu_0", fitFunction = fitFn, predFunction = predFn, currentFit = list(mu1 = 0, mu2 = 1, KL_div = 1), family = 'binomial'),
  learner3 = VEBBoostNode$new("mu_0", fitFunction = fitFn, predFunction = predFn, currentFit = list(mu1 = 0, mu2 = 1, KL_div = 1), family = 'binomial'),
  learner4 = VEBBoostNode$new("mu_0", fitFunction = fitFn, predFunction = predFn, currentFit = list(mu1 = 0, mu2 = 1, KL_div = 1), family = 'binomial'),
  learner5 = VEBBoostNode$new("mu_0", fitFunction = fitFn, predFunction = predFn, currentFit = list(mu1 = 0, mu2 = 1, KL_div = 1), family = 'binomial'),
  learner6 = VEBBoostNode$new("mu_0", fitFunction = fitFn, predFunction = predFn, currentFit = list(mu1 = 0, mu2 = 1, KL_div = 1), family = 'binomial'),
  learner7 = VEBBoostNode$new("mu_0", fitFunction = fitFn, predFunction = predFn, currentFit = list(mu1 = 0, mu2 = 1, KL_div = 1), family = 'binomial'),
  learner8 = VEBBoostNode$new("mu_0", fitFunction = fitFn, predFunction = predFn, currentFit = list(mu1 = 0, mu2 = 1, KL_div = 1), family = 'binomial'),
  learner9 = VEBBoostNode$new("mu_0", fitFunction = fitFn, predFunction = predFn, currentFit = list(mu1 = 0, mu2 = 1, KL_div = 1), family = 'binomial')
)

for (learner in learner_MNIST$learners) {
  learner$ensemble = learner_MNIST
}

learner_MNIST$learners$learner0$Y = as.numeric(Y == 0)
learner_MNIST$learners$learner1$Y = as.numeric(Y == 1)
learner_MNIST$learners$learner2$Y = as.numeric(Y == 2)
learner_MNIST$learners$learner3$Y = as.numeric(Y == 3)
learner_MNIST$learners$learner4$Y = as.numeric(Y == 4)
learner_MNIST$learners$learner5$Y = as.numeric(Y == 5)
learner_MNIST$learners$learner6$Y = as.numeric(Y == 6)
learner_MNIST$learners$learner7$Y = as.numeric(Y == 7)
learner_MNIST$learners$learner8$Y = as.numeric(Y == 8)
learner_MNIST$learners$learner9$Y = as.numeric(Y == 9)

tol = 4.2

system.time({
  learner_MNIST$convergeFit(tol)
})

system.time({
  learner_MNIST$convergeFitAll(tol)
})

system.time({
  while ((abs(tail(tail(learner_MNIST$ELBO_progress, 1)[[1]], 1) - tail(tail(learner_MNIST$ELBO_progress, 2)[[1]], 1)) > tol) && (sum(sapply(learner_MNIST$learners, function(x) length(Traverse(x, filterFun = function(node) node$isLeaf & !node$isLocked)))) > 0)) {
    learner_MNIST$convergeFitAll(tol)
  }
})


test <- read.csv("D:/School/Fall 2019/Matthew Stephens Research/Kaggle MNIST Data/test.csv", stringsAsFactors=FALSE)
X_new = as.matrix(test)
X_new = X_new / 255
X_new_stumps = make_stumps_matrix(X_new, F, do.call(cbind, sapply(X_stumps, function(x) attr(x, 'br'))))


learner_MNIST$predict.veb(X_new_stumps, 1)

preds = apply(learner_MNIST$pred_mu1, MARGIN = 1, which.max) - 1

submission = data.frame(ImageId = 1:length(preds), Label = preds)

write.csv(submission, "D:/School/Fall 2019/Matthew Stephens Research/Kaggle MNIST Data/submission.csv", row.names = F)
