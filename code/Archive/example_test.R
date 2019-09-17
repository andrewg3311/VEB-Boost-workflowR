set.seed(1138)

X1 = rnorm(100)
X2 = rnorm(100)
X3 = rnorm(100)
X4 = rnorm(100)
X5 = rnorm(100)

beta1 = rnorm(1)
beta2 = rnorm(1)
beta3 = rnorm(1)
beta4 = rnorm(1)
beta5 = rnorm(1)

Y = ((X1*beta1) * ((X2*beta2) + X3*beta3)) + ((X4*beta4) * X5*beta5) + rnorm(100)

fitFn1 = function(Y, sigma2, init) {
  fit = lm(Y ~ X1, weights = 1/sqrt(sigma2))
  mu1 = predict(fit)
  mu2 = mu1^2 + (X1^2 * summary(fit)$coef[2, "Std. Error"]^2)
  KL_div = .1
  
  return(list(fit = fit, mu1 = mu1, mu2 = mu2, KL_div = KL_div))
}

fitFn2 = function(Y, sigma2, init) {
  fit = lm(Y ~ X2, weights = 1/sqrt(sigma2))
  mu1 = predict(fit)
  mu2 = mu1^2 + (X2^2 * summary(fit)$coef[2, "Std. Error"]^2)
  KL_div = .1
  
  return(list(fit = fit, mu1 = mu1, mu2 = mu2, KL_div = KL_div))
}

fitFn3 = function(Y, sigma2, init) {
  fit = lm(Y ~ X3, weights = 1/sqrt(sigma2))
  mu1 = predict(fit)
  mu2 = mu1^2 + (X3^2 * summary(fit)$coef[2, "Std. Error"]^2)
  KL_div = .1
  
  return(list(fit = fit, mu1 = mu1, mu2 = mu2, KL_div = KL_div))
}

fitFn4 = function(Y, sigma2, init) {
  fit = lm(Y ~ X4, weights = 1/sqrt(sigma2))
  mu1 = predict(fit)
  mu2 = mu1^2 + (X4^2 * summary(fit)$coef[2, "Std. Error"]^2)
  KL_div = .1
  
  return(list(fit = fit, mu1 = mu1, mu2 = mu2, KL_div = KL_div))
}

fitFn5 = function(Y, sigma2, init) {
  fit = lm(Y ~ X5, weights = 1/sqrt(sigma2))
  mu1 = predict(fit)
  mu2 = mu1^2 + (X5^2 * summary(fit)$coef[2, "Std. Error"]^2)
  KL_div = .1
  
  return(list(fit = fit, mu1 = mu1, mu2 = mu2, KL_div = KL_div))
}


learner = VEBBoostNode$new(".", operator = "+")
L = learner$AddChildVEB("L", operator = "*")$
  AddChildVEB("mu_1", fitFunction = fitFn1, currentFit = fitFn1(Y, rep(1, 100), NULL))$parent$
  AddChildVEB("LR", operator = "+")$
  AddChildVEB("mu_2", fitFunction = fitFn2, currentFit = fitFn2(Y, rep(1, 100), NULL))$parent$
  AddChildVEB("mu_3", fitFunction = fitFn3, currentFit = fitFn3(Y, rep(1, 100), NULL))$root$
  AddChildVEB("R", operator = "*")$
  AddChildVEB("mu_4", fitFunction = fitFn4, currentFit = fitFn4(Y, rep(1, 100), NULL))$parent$
  AddChildVEB("mu_5", fitFunction = fitFn5, currentFit = fitFn5(Y, rep(1, 100), NULL))

learner$Y = Y
learner$sigma2 = 1


for (j in 1:10) {learner$Do(function(x) x$updateFit(), filterFun = function(x) x$isLeaf)}