
# predFn2_stumps = function(X_new, currentFit, moment = c(1, 2)) {
#   X_avg_mat_new = matrix(currentFit$X_avg, nrow = nrow(X_new), ncol = ncol(X_new), byrow = T)
#   X_cent_new = X_new - X_avg_mat_new
#   beta_post_1 = currentFit$alpha * currentFit$mu
#   if (moment == 1) {
#     return(as.numeric(currentFit$Y_avg + (X_cent_new %*% (beta_post_1))))
#   } else if (moment == 2) {
#     beta_post_2 = currentFit$alpha * (currentFit$sigma2 + currentFit$mu^2)
#     return(as.numeric((currentFit$Y_avg^2) + 2*currentFit$Y_avg*(X_cent_new %*% beta_post_1) + ((X_cent_new^2) %*% beta_post_2)))
#   } else {
#     stop("`moment` must be either 1 or 2")
#   }
# }

# could still even be more efficient, but what ever
predFn2_stumps = function(X_new, currentFit, moment = c(1, 2)) {
  beta_post_1 = currentFit$alpha * currentFit$mu
  beta_post_1_mat = matrix(beta_post_1, ncol = ncol(X_new), byrow = F)
  beta_post_1_rev_mat = apply(beta_post_1_mat, MARGIN = 2, rev)
  X = attr(X_new, 'X') # original X stored in X_new
  if (moment == 1) {
    pred_mu1 = rep(currentFit$intercept, nrow(X_new))
    for (i in 1:nrow(X_new)) {
      for (j in 1:ncol(X)) {
        pred_mu1[i] = pred_mu1[i] + sum((X_new[i, j] >= X[attr(X, 'order')[, j], j]) * beta_post_1_rev_mat[, j])
      }
    }
    return(pred_mu1)
  } else if (moment == 2) {
    stop("moment = 2 not yet implemented")
  } else {
    stop("`moment` must be either 1 or 2")
  }
}


# let attr(X, 'order') be an nxp matrix, each column is the DECREASING order of variable X_j

weighted_SER2_stumps = function(X, Y, sigma2, init = list(V = NULL)) {
  if (length(sigma2) == 1) {
    sigma2 = rep(sigma2, nrow(X))
  }
  
  w = (1 / sigma2) / sum(1 / sigma2)
  #X_avg = rev(cumsum(w[order(X, decreasing = T)]))
  X_avg = apply(attr(X, 'order'), MARGIN = 2, function(ord) rev(cumsum(w[ord]))) # each column of X_avg is column mean for that var's stump matrix
  #Y_avg = weighted.mean(Y, 1/sigma2)
  Y_avg = sum(Y * w)
  Y_cent = Y - Y_avg
  
  prior_weights = rep(1 / prod(dim(X)), prod(dim(X)))
  #tau_no_V = (t(X_cent^2) %*% (1 / sigma2))
  tau_no_V = as.numeric(sapply(1:ncol(attr(X, 'order')), function(j) rev(cumsum((1 / sigma2)[attr(X, 'order')[, j]])) - 2*rev(cumsum((1 / sigma2)[attr(X, 'order')[, j]]))*X_avg[, j] + (X_avg[, j]^2 * sum(1 / sigma2))))
  #nu = colSums(X_cent * Y_cent / sigma2)
  nu = as.numeric(sapply(1:ncol(attr(X, 'order')), function(j) rev(cumsum((Y_cent / sigma2)[attr(X, 'order')[, j]])) - (X_avg[, j] * sum(Y_cent / sigma2))))
  
  V = ifelse(is.null(init$V), 1, init$V)
  V = optimize_V(tau_no_V, nu, sigma2, prior_weights, V)
  
  tau = tau_no_V + (1 / V)
  
  alpha = log(prior_weights) - (.5 * log(tau)) + (.5 * nu^2 / tau)
  alpha = alpha - max(alpha)
  alpha = exp(alpha)
  alpha = alpha / sum(alpha)
  
  mu = nu / tau
  
  sigma2_post = 1 / tau
  
  beta_post_1 = alpha * mu
  beta_post_2 = alpha * (sigma2_post + mu^2)
  intercept = as.numeric(Y_avg - as.numeric(X_avg)%*%beta_post_1)
  
  # mu1 = predFn2_stumps(X, list(alpha = alpha, mu = mu, sigma2_post = sigma2_post, intercept = intercept, X_avg = X_avg, Y_avg = Y_avg), 1)
  # mu2 = predFn2_stumps(X, list(alpha = alpha, mu = mu, sigma2_post = sigma2_post, intercept = intercept, X_avg = X_avg, Y_avg = Y_avg), 2)
  
  # mu1 = predFn2(attr(X, 'X_stumps'), list(alpha = matrix(alpha, ncol = 1), mu = matrix(mu, ncol = 1), sigma2_post = matrix(sigma2_post, ncol = 1), intercept = intercept, X_avg = as.numeric(X_avg), Y_avg = Y_avg), 1)
  # mu2 = predFn2(attr(X, 'X_stumps'), list(alpha = matrix(alpha, ncol = 1), mu = matrix(mu, ncol = 1), sigma2_post = matrix(sigma2_post, ncol = 1), intercept = intercept, X_avg = as.numeric(X_avg), Y_avg = Y_avg), 2)
  
  beta_1_mat = matrix(beta_post_1, ncol = ncol(X), byrow = F)
  beta_2_mat = matrix(beta_post_2, ncol = ncol(X), byrow = F)
  mu1 = intercept + rowSums(sapply(1:ncol(X), function(j) cumsum(beta_1_mat[, j])[attr(X, 'rank')[, j]]))
  mu2 = Y_avg^2 + 2*Y_avg*(rowSums(sapply(1:ncol(X), function(j) cumsum(beta_1_mat[, j])[attr(X, 'rank')[, j]])) - sum(X_avg*beta_1_mat)) +
    rowSums(sapply(1:ncol(X), function(j) cumsum(beta_2_mat[, j])[attr(X, 'rank')[, j]])) - 2*rowSums(sapply(1:ncol(X), function(j) cumsum(beta_2_mat[, j] * X_avg[, j])[attr(X, 'rank')[, j]])) + sum(X_avg^2 * beta_2_mat)
  
  # if (any(mu2 < mu1^2)) { # sanity check, this should be >= 0
  #   browser()
  #   stop("Predicted variance is negative")
  # }
  
  KL_div = calc_KL(matrix(mu, ncol = 1), matrix(alpha, ncol = 1), matrix(sigma2_post, ncol = 1), V)
  
  return(list(mu1 = as.numeric(mu1), mu2 = as.numeric(mu2), KL_div = KL_div, alpha = alpha, mu = mu, sigma2_post = sigma2_post, intercept = intercept, V = V, X_avg = X_avg, Y_avg = Y_avg))
}

fitFn2_stumps = function(X, Y, sigma2, init) {
  res = weighted_SER2_stumps(X, Y, sigma2, init)
  return(res)
}

# this version includes a distribution for the intercept, N(0, s_i^2), fit as independent, do EB for s_i^2
weighted_SER2_stumps_int = function(X, Y, sigma2, init = list(V = NULL, intercept = NULL)) {
  if (length(sigma2) == 1) {
    sigma2 = rep(sigma2, nrow(X))
  }
  
  sum_inv_sigma2 = sum(1 / sigma2)
  
  intercept = ifelse(is.null(init$intercept), 0, init$intercept)
  Y_cent = Y - intercept
  
  prior_weights = rep(1 / prod(dim(X)), prod(dim(X)))
  tau_no_V = as.numeric(sapply(1:ncol(attr(X, 'order')), function(j) rev(cumsum((1 / sigma2)[attr(X, 'order')[, j]]))))
  nu = as.numeric(sapply(1:ncol(attr(X, 'order')), function(j) rev(cumsum((Y_cent / sigma2)[attr(X, 'order')[, j]]))))
  
  V = ifelse(is.null(init$V), 1, init$V)
  V = optimize_V(tau_no_V, nu, sigma2, prior_weights, V)
  
  tau = tau_no_V + (1 / V)
  
  alpha = log(prior_weights) - (.5 * log(tau)) + (.5 * nu^2 / tau)
  alpha = alpha - max(alpha)
  alpha = exp(alpha)
  alpha = alpha / sum(alpha)
  
  mu = nu / tau
  
  sigma2_post = 1 / tau
  
  beta_post_1 = alpha * mu
  beta_post_2 = alpha * (sigma2_post + mu^2)
  
  beta_1_mat = matrix(beta_post_1, ncol = ncol(X), byrow = F)
  beta_2_mat = matrix(beta_post_2, ncol = ncol(X), byrow = F)
  mu1 = rowSums(sapply(1:ncol(X), function(j) cumsum(beta_1_mat[, j])[attr(X, 'rank')[, j]]))
  mu2 = rowSums(sapply(1:ncol(X), function(j) cumsum(beta_2_mat[, j])[attr(X, 'rank')[, j]]))
  
  # EB for intercept
  Y_int = Y - mu1
  nu_int = sum(Y_int / sigma2)
  V_int = max(0, ((nu_int^2 / sum_inv_sigma2) - 1) / sum_inv_sigma2)
  tau_int = (1 / V_int) + sum_inv_sigma2
  intercept = ifelse(V_int == 0, 0, nu_int / tau_int)
  mu2_int = ifelse(V_int == 0, 0, 1 / tau_int + intercept^2)
  

  # if (any(mu2 < mu1^2)) { # sanity check, this should be >= 0
  #   browser()
  #   stop("Predicted variance is negative")
  # }
  
  mu2 = mu2 + 2*mu1*intercept + mu2_int
  mu1 = mu1 + intercept
  
  KL_div_int = ifelse(V_int == 0, 0, ((log(V_int) + log(tau_int)) + (((1 / tau_int) + (intercept^2)) /  V_int) - 1) / 2) # KL div from Q_int to prior
  
  KL_div = calc_KL(matrix(mu, ncol = 1), matrix(alpha, ncol = 1), matrix(sigma2_post, ncol = 1), V) + KL_div_int
  
  return(list(mu1 = as.numeric(mu1), mu2 = as.numeric(mu2), KL_div = KL_div, alpha = alpha, mu = mu, sigma2_post = sigma2_post, intercept = intercept, V = V, V_int = V_int))
}

fitFn2_stumps_int = function(X, Y, sigma2, init) {
  res = weighted_SER2_stumps_int(X, Y, sigma2, init)
  return(res)
}


make_X_new_stumps = function(X, X_new) {
  X_new_stumps = matrix(0, nrow = nrow(X_new), ncol = prod(dim(X)))
  for (i in 1:nrow(X_new)) {
    row = numeric()
    for (j in 1:ncol(X)) {
      row = c(row, 1 * (X_new[i, j] >= rev(X[attr(X, 'order')[, j], j])))
    }
    X_new_stumps[i, ] = row
  }
  return(X_new_stumps)
}


weighted_SER2_lin_stumps = function(X, Y, sigma2, init = list(V = NULL)) {
  if (length(sigma2) == 1) {
    sigma2 = rep(sigma2, nrow(X))
  }
  
  w = (1 / sigma2) / sum(1 / sigma2)
  #X_avg = rev(cumsum(w[order(X, decreasing = T)]))
  X_avg_lin = apply(X, MARGIN = 2, function(col) weighted.mean(col, 1/sigma2))
  X_avg_stumps = apply(attr(X, 'order'), MARGIN = 2, function(ord) rev(cumsum(w[ord]))) # each column of X_avg is column mean for that var's stump matrix
  #Y_avg = weighted.mean(Y, 1/sigma2)
  Y_avg = sum(Y * w)
  Y_cent = Y - Y_avg
  X_avg_mat_lin = matrix(X_avg_lin, nrow = nrow(X), ncol = ncol(X), byrow = T)
  X_cent_lin = X - X_avg_mat_lin
  
  prior_weights = rep(1 / (ncol(X) + prod(dim(X))), (ncol(X) + prod(dim(X))))
  #tau_no_V = (t(X_cent^2) %*% (1 / sigma2))
  tau_no_V_lin = (t(X_cent_lin^2) %*% (1 / sigma2))
  nu_lin = colSums(X_cent_lin * Y / sigma2)
  tau_no_V_stumps = as.numeric(sapply(1:ncol(attr(X, 'order')), function(j) rev(cumsum((1 / sigma2)[attr(X, 'order')[, j]])) - 2*rev(cumsum((1 / sigma2)[attr(X, 'order')[, j]]))*X_avg_stumps[, j] + (X_avg_stumps[, j]^2 * sum(1 / sigma2))))
  #nu = colSums(X_cent * Y_cent / sigma2)
  nu_stumps = as.numeric(sapply(1:ncol(attr(X, 'order')), function(j) rev(cumsum((Y_cent / sigma2)[attr(X, 'order')[, j]])) - (X_avg_stumps[, j] * sum(Y_cent / sigma2))))
  
  tau_no_V = c(tau_no_V_lin, tau_no_V_stumps)
  nu = c(nu_lin, nu_stumps)
  
  V = ifelse(is.null(init$V), 1, init$V)
  V = optimize_V(tau_no_V, nu, sigma2, prior_weights, V)
  
  tau = tau_no_V + (1 / V)
  
  alpha = log(prior_weights) - (.5 * log(tau)) + (.5 * nu^2 / tau)
  alpha = alpha - max(alpha)
  alpha = exp(alpha)
  alpha = alpha / sum(alpha)
  
  mu = nu / tau
  
  sigma2_post = 1 / tau
  
  beta_post_1 = alpha * mu
  beta_post_2 = alpha * (sigma2_post + mu^2)
  beta_post_1_lin = beta_post_1[1:ncol(X)]
  beta_post_1_stumps = beta_post_1[-c(1:ncol(X))]
  beta_post_2_lin = beta_post_2[1:ncol(X)]
  beta_post_2_stumps = beta_post_2[-c(1:ncol(X))]
  intercept = as.numeric(Y_avg - c(X_avg_lin, as.numeric(X_avg_stumps))%*%beta_post_1)
  
  # mu1 = predFn2_stumps(X, list(alpha = alpha, mu = mu, sigma2_post = sigma2_post, intercept = intercept, X_avg = X_avg, Y_avg = Y_avg), 1)
  # mu2 = predFn2_stumps(X, list(alpha = alpha, mu = mu, sigma2_post = sigma2_post, intercept = intercept, X_avg = X_avg, Y_avg = Y_avg), 2)
  
  # mu1 = predFn2(attr(X, 'X_stumps'), list(alpha = matrix(alpha, ncol = 1), mu = matrix(mu, ncol = 1), sigma2_post = matrix(sigma2_post, ncol = 1), intercept = intercept, X_avg = as.numeric(X_avg), Y_avg = Y_avg), 1)
  # mu2 = predFn2(attr(X, 'X_stumps'), list(alpha = matrix(alpha, ncol = 1), mu = matrix(mu, ncol = 1), sigma2_post = matrix(sigma2_post, ncol = 1), intercept = intercept, X_avg = as.numeric(X_avg), Y_avg = Y_avg), 2)
  
  beta_1_mat_stumps = matrix(beta_post_1_stumps, ncol = ncol(X), byrow = F)
  beta_2_mat_stumps = matrix(beta_post_2_stumps, ncol = ncol(X), byrow = F)
  mu1 = intercept + as.numeric(X %*% beta_post_1_lin) + rowSums(sapply(1:ncol(X), function(j) cumsum(beta_1_mat_stumps[, j])[attr(X, 'rank')[, j]]))
  mu2 = Y_avg^2 + 2*Y_avg*((X_cent_lin %*% beta_post_1_lin) + rowSums(sapply(1:ncol(X), function(j) cumsum(beta_1_mat_stumps[, j])[attr(X, 'rank')[, j]])) - sum(X_avg_stumps*beta_1_mat_stumps)) + ((X_cent_lin^2) %*% beta_post_2_lin) + 
    rowSums(sapply(1:ncol(X), function(j) cumsum(beta_2_mat_stumps[, j])[attr(X, 'rank')[, j]])) - 2*rowSums(sapply(1:ncol(X), function(j) cumsum(beta_2_mat_stumps[, j] * X_avg_stumps[, j])[attr(X, 'rank')[, j]])) + sum(X_avg_stumps^2 * beta_2_mat_stumps)
  
  # if (any(mu2 < mu1^2)) { # sanity check, this should be >= 0
  #   browser()
  #   stop("Predicted variance is negative")
  # }
  
  KL_div = calc_KL(matrix(mu, ncol = 1), matrix(alpha, ncol = 1), matrix(sigma2_post, ncol = 1), V)
  
  return(list(mu1 = as.numeric(mu1), mu2 = as.numeric(mu2), KL_div = KL_div, alpha = alpha, mu = mu, sigma2_post = sigma2_post, intercept = intercept, V = V, X_avg = c(X_avg_lin, as.numeric(X_avg_stumps)), Y_avg = Y_avg))
}


# could definitely be faster
predFn_lin_stumps = function(X_new, currentFit, moment = c(1, 2)) {
    beta_post_1 = currentFit$alpha * currentFit$mu
    beta_post_1_lin = beta_post_1[1:ncol(X_new)]
    beta_post_1_stumps = beta_post_1[-c(1:ncol(X_new))]
    beta_post_1_stumps_mat = matrix(beta_post_1_stumps, ncol = ncol(X_new), byrow = F)
    beta_post_1_stumps_cumsum_mat = apply(beta_post_1_stumps_mat, MARGIN = 2, cumsum)
    X = attr(X_new, 'X') # original X stored in X_new
    X_sorted = apply(sapply(1:ncol(X), function(j) X[attr(X, 'order')[, j], j]), MARGIN = 2, rev)
    if (moment == 1) {
      pred_mu1 = rep(currentFit$intercept, nrow(X_new))
      for (i in 1:nrow(X_new)) {
        #pred_mu1[i] = pred_mu1[i] + sum((matrix(X_new[i, ], nrow = nrow(X_sorted), ncol = ncol(X_sorted), byrow = T) >= X_sorted) * beta_post_1_stumps_mat)
        positions = sapply(1:ncol(X_new), function(j) findInterval(X_new[i, j], X_sorted[, j])) # findInterval uses binary search in sorted list
        pred_mu1[i] = pred_mu1[i] + sum(diag(as.matrix(beta_post_1_stumps_cumsum_mat[positions[which(positions != 0)], which(positions != 0)])))
      }
      pred_mu1 = pred_mu1 + (X_new %*% beta_post_1_lin)
      return(as.numeric(pred_mu1))
  } else if (moment == 2) {
    stop("moment = 2 not yet implemented")
  } else {
    stop("`moment` must be either 1 or 2")
  }
}



fitFn2_lin_stumps = function(X, Y, sigma2, init) {
  res = weighted_SER2_lin_stumps(X, Y, sigma2, init)
  return(res)
}