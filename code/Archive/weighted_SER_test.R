set.seed(1138)

n = 1000

X1 = rnorm(n)
X2 = rnorm(n)
X3 = rnorm(n)
X4 = rnorm(n)
X5 = rnorm(n)
X1_null = rnorm(n)
X2_null = rnorm(n)
X3_null = rnorm(n)
X4_null = rnorm(n)
X5_null = rnorm(n)

beta1 = rnorm(1)
beta2 = rnorm(1)
beta3 = rnorm(1)
beta4 = rnorm(1)
beta5 = rnorm(1)

mu_true = ((X1*beta1) * ((X2*beta2) + (X3*beta3))) + ((X4*beta4) * (X5*beta5)) 

Y = mu_true + rnorm(n)

calc_Q = function(X, Sigma2, Mu, Alpha) {
  # X is data matrix
  # Sigma2[j, l] is the posterior variance for b_l when entry j is selected, p x L
  # Mu[j, l] is the posterior mean for b_l when entry j is selected, p x L
  # Alpha[j, l] is the posterior probability selecting entry j from b_l, p x L
  # Z is matrix of covariates (e.g. column for intercept, top 10 PCs, etc)
  # delta is current estimate for effects of Z variables
  
  ASU2 = Alpha * (Sigma2 + Mu^2) # [j, l] = alpha[j, l] * (Sigma2[j, l] + Mu[j, l]^2)
  
  Q = rowSums(X^2 %*% ASU2) # start w/ sum_l E_(ql)[(x_i' b_l)^2]
  
  return(Q)
  
}


calc_KL = function(Mu, Alpha, Sigma2, prior_var = 1) {
  p = nrow(Mu)
  L = ncol(Mu)
  prior_weights = rep(1/p, p)
  P = matrix(prior_weights, nrow = p, ncol = L)
  b_post = rowSums(Alpha * Mu)
  
  prior_var = matrix(prior_var, nrow = p, ncol = L, byrow = T)
  
  KL_div = Alpha * (log(Alpha) - log(P) + (log(prior_var) / 2) - (log(Sigma2) / 2) - .5 + ((Sigma2 + Mu^2) / (2 * prior_var)))
  KL_div[Alpha == 0] = 0
  return(sum(KL_div))
}

weighted_SER = function(X, Y, sigma2, prior_var = 1) {
  if (length(sigma2) == 1) {
    sigma2 = rep(sigma2, nrow(X))
  }
  prior_weights = rep(1 / ncol(X), ncol(X))
  tau = (t(X^2) %*% (1 / sigma2)) + (1 / prior_var)
  nu = colSums(X * Y / sigma2)
  
  alpha = log(prior_weights) - (.5 * log(tau)) + (.5 * nu^2 / tau)
  alpha = alpha - max(alpha)
  alpha = exp(alpha)
  alpha = alpha / sum(alpha)
  
  mu = nu / tau
  
  sigma2_post = 1 / tau
  
  mu1 = X %*% (alpha * mu)
  mu2 = calc_Q(X, sigma2_post, mu, alpha)
  KL_div = calc_KL(mu, alpha, sigma2_post, prior_var)
  
  return(list(mu1 = as.numeric(mu1), mu2 = as.numeric(mu2), KL_div = KL_div, alpha = alpha, mu = mu, sigma2_post = sigma2_post))
}

fitFn1 = function(Y, sigma2, init) {
  return(weighted_SER(X = cbind(X1, X1_null), Y, sigma2))
}

fitFn2 = function(Y, sigma2, init) {
  return(weighted_SER(X = cbind(X2, X2_null), Y, sigma2))
}

fitFn3 = function(Y, sigma2, init) {
  return(weighted_SER(X = cbind(X3, X3_null), Y, sigma2))
}

fitFn4 = function(Y, sigma2, init) {
  return(weighted_SER(X = cbind(X4, X4_null), Y, sigma2))
}

fitFn5 = function(Y, sigma2, init) {
  return(weighted_SER(X = cbind(X5, X5_null), Y, sigma2))
}


learner = VEBBoostNode$new(".", operator = "+")
learner$Y = Y
learner$sigma2 = 1
L = learner$AddChildVEB("L", operator = "*")$
  AddChildVEB("mu_1", fitFunction = fitFn1, currentFit = fitFn1(Y, rep(1, n), NULL))$parent$
  AddChildVEB("LR", operator = "+")$
  AddChildVEB("mu_2", fitFunction = fitFn2, currentFit = fitFn2(Y, rep(1, n), NULL))$parent$
  AddChildVEB("mu_3", fitFunction = fitFn3, currentFit = fitFn3(Y, rep(1, n), NULL))$root$
  AddChildVEB("R", operator = "*")$
  AddChildVEB("mu_4", fitFunction = fitFn4, currentFit = fitFn4(Y, rep(1, n), NULL))$parent$
  AddChildVEB("mu_5", fitFunction = fitFn5, currentFit = fitFn5(Y, rep(1, n), NULL))

ELBOs = numeric(100)
ELBOs[1] = -Inf
ELBOs[2] = learner$ELBO
i = 2
while (ELBOs[i] - ELBOs[i-1] > 1e-6) {
  learner$Do(function(x) x$updateFit(), filterFun = function(x) x$isLeaf)
  i = i+1
  ELBOs[i] = learner$ELBO
}
ELBOs = ELBOs[2:i]


# 
# learner2 = VEBBoostNode$new(".", operator = "+")
# learner2$Y = Y
# learner2$sigma2 = 1
# L2 = learner2$AddChildVEB("L", operator = "*")$
#   AddChildVEB("mu_1", fitFunction = fitFn1, currentFit = list(mu1 = rep(0, 100), mu2 = rep(1, 100), KL_div = 0))$parent$
#   AddChildVEB("LR", operator = "+")$
#   AddChildVEB("mu_2", fitFunction = fitFn2, currentFit = list(mu1 = rep(0, 100), mu2 = rep(1, 100), KL_div = 0))$parent$
#   AddChildVEB("mu_3", fitFunction = fitFn3, currentFit = list(mu1 = rep(0, 100), mu2 = rep(1, 100), KL_div = 0))$root$
#   AddChildVEB("R", operator = "*")$
#   AddChildVEB("mu_4", fitFunction = fitFn4, currentFit = list(mu1 = rep(0, 100), mu2 = rep(1, 100), KL_div = 0))$parent$
#   AddChildVEB("mu_5", fitFunction = fitFn5, currentFit = list(mu1 = rep(0, 100), mu2 = rep(1, 100), KL_div = 0))
# 
# 
# ELBOs2 = numeric(100)
# ELBOs2[1] = -Inf
# ELBOs2[2] = learner2$ELBO
# i = 2
# while (ELBOs2[i] - ELBOs2[i-1] > 1e-6) {
#   learner2$Do(function(x) x$updateFit(), filterFun = function(x) x$isLeaf)
#   i = i+1
#   ELBOs2[i] = learner2$ELBO
# }
# ELBOs2 = ELBOs2[2:i]