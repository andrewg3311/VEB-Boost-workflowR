library(data.tree)
# get SuSiE stumps source code
sapply(list.files("D:/School/Fall 2019/Matthew Stephens Research/susieR-susie_stumps/R/"), function(x) source(paste("D:/School/Fall 2019/Matthew Stephens Research/susieR-susie_stumps/R/", x, sep = "")))

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

# over-write optimize_V code w/ weighted version
log_lik_SER = function(V, tau_no_V, nu, sigma2, prior_weights) {
  if (length(sigma2) == 1) {
    sigma2 = rep(sigma2, length(nu))
  }
  tau = tau_no_V + (1 / V)
  m = -(log(tau) / 2) + (nu^2 / (2 * tau))
  m_max = max(m)
  w = exp(m - m_max)
  -(log(V) / 2) + m_max + log(sum(prior_weights * w))
}

neg.loglik.logscale = function(lV, tau_no_V, nu, sigma2, prior_weights){
  -log_lik_SER(exp(lV), tau_no_V, nu, sigma2, prior_weights)
}

optimize_V = function(tau_no_V, nu, sigma2, prior_weights, V = 1) {
  lV = optim(par = log(V), fn = neg.loglik.logscale, tau_no_V = tau_no_V, nu = nu, sigma2 = sigma2, prior_weights = prior_weights, method='Brent', lower = -10, upper = 15)$par
  V = exp(lV)
  return(V)
}



library(doParallel)
### MAKE PARELLE VERSIONS OF MATRIX MULTIPLICATION ###
#' @title Computes standardized.X \%*\% b
#' @param X an n by p matrix with three attributes: scaled:center, scaled:scale, and attr(X, 'd')
#' @param b a p vector
#' @return an n vector
#' @importFrom Matrix t
#' @importFrom Matrix tcrossprod
#' @keywords internal
compute_Xb_par = function(X, b){
  if(is.list(X)){
    n_var = unlist(lapply(X,get_ncol)) # number of variables for each element of X
    b_split = split_vector(b,n_var) # split b into a list of vectors
    #Xb = mapply(compute_Xb,X,b_split,SIMPLIFY=FALSE) # apply compute_Xb to elements of lists X,b_split
    Xb = foreach(i = 1:length(X), .combine = '+', .inorder = F) %dopar% {
      compute_Xb(X[[i]], b_split[[i]])
    }
    return(Xb) # add the results up
  } else {
    cm = get_cm(X)
    csd = get_csd(X)
    #scale Xb
    #when X is a trend filtering matrix or stumps matrix, use special matrix mult
    if (is.tfmatrix(X))
      scaled.Xb <- compute_tf_Xb(get_order(X), b/csd)
    else if(is.stumps_matrix(X)){
      scaled.Xb <- compute_stumps_Xb(attr(X, 'Xord'),  b/csd)
    } else if(is.tfg_matrix(X)){
      scaled.Xb <- compute_tfg_Xb(X,b/csd)
    } else
      #when X is an ordinary sparse/dense matrix
      scaled.Xb <- tcrossprod(X, t(b/csd))
    #center Xb
    Xb <- scaled.Xb - sum(cm*b/csd)
    return(as.numeric(Xb))
  }
}

#' @title Computes t(standardized.X)\%*\%y using sparse multiplication trick
#' @param X an n by p unstandardized matrix with three attributes: scaled:center, scaled:scale, and attr(X, 'd')
#' @param y an n vector
#' @return a p vector
#' @importFrom Matrix t
#' @importFrom Matrix crossprod
compute_Xty_par = function(X, y){
  if(is.list(X)){ # perform Xty for each element in list and concatenate the results
    # unlist(lapply(X,compute_Xty,y=y))
    return(
      foreach(Xj = X, .combine = 'c') %dopar% {
        compute_Xty(Xj, y)
      }
    )
  } else {
    cm = get_cm(X)
    csd = get_csd(X)

    #when X is a trend filtering matrix
    if (is.tfmatrix(X))
      scaled.Xty <- compute_tf_Xty(get_order(X),y)/csd
    else if(is.stumps_matrix(X))
      scaled.Xty <- compute_stumps_Xty(attr(X,'Xord'),y)/csd
    else if(is.tfg_matrix(X))
      scaled.Xty <- compute_tfg_Xty(X,y)/csd
    #when X is an ordinary sparse/dense matrix
    else{
      ytX <- crossprod(y, X)
      scaled.Xty <- t(ytX/csd)
    }
    #center Xty
    centered.scaled.Xty <- scaled.Xty - cm/csd * sum(y)
    return(as.numeric(centered.scaled.Xty))
  }
}




# over-write make_stumps_matrix to not rely on susieR::
# also change Xtrain to be a list (allows for different lengths of breaks)
make_stumps_matrix = function(X, include_linear, Xtrain=NULL){
  if(is.null(Xtrain)){Xtrain = lapply(1:ncol(X), function(i) X[, i])}
  
  xl=list() # initialize
  if(include_linear){ #include X as a regular matrix first
    attr(X,"nrow") <- nrow(X)
    attr(X,"ncol") <- ncol(X)
    attr(X,"scaled:center") <- rep(0,ncol(X))
    attr(X,"scaled:scale") <- rep(1,ncol(X))
    xl=c(xl,list(X))
  }
  
  for(i in 1:ncol(X)){xl= c(xl,list(make_tfg_matrix(X[,i],Xtrain[[i]])))}
  return(xl)
}

# over-write stumps multiplication (WHY IS THERE A -1*.... ?!?!?)
#' @title Compute unscaled X \%*\% b using the special structure of trend filtering
#' @param X a tfg_matrix created by make_tfg_matrix
#' @param b a p vector of the changes at each change point
#' @return an n vector of the means at each data point
#' @keywords internal
compute_tfg_Xb = function(X,b){
  order = get_order(X)
  for(i in 1:(order+1)){
    #b = rev(cumsum(rev(b))) # computes mean in each bin
    b = spatstat.utils::revcumsum(b) # faster than rev(cumsum(rev(b)))
  }
  return(b[attr(X,"t_to_bin")]) #  maps bin means to a mean for each datapoint
}

#' @title Compute t(X) \%*\% y using the special structure of trend filtering
#' @param X a tfg_matrix created by make_tfg_matrix
#' @param y an n vector of data
#' @return a p vector
#' @keywords internal
compute_tfg_Xty = function(X,y){
  order = get_order(X)
  y = y[attr(X,"order_t")] # sort y according to increasing t
  for (i in 1:(order+1)){
    y = cumsum(y)
  }
  return(y[attr(X,"bin_to_t")])
}

# computes (X - X_avg)^2 %*% b
compute_X2b = function(X, b, X_avg = 0) {
  if (is.list(X)) {
    n_var = unlist(lapply(X,get_ncol)) # number of variables for each element of X
    b_split = split_vector(b,n_var) # split b into a list of vectors
    X2b = mapply(compute_X2b, X, b_split, SIMPLIFY = FALSE) # apply compute_X2b to elements of lists X, b_split
    return(Reduce(`+`, X2b) - 2*compute_Xb(X, b*X_avg) + sum(X_avg^2 * b)) # add the results up
  } else {
    if (is.tfg_matrix(X)) {
      # X is boolean matrix, so X^2 = X
      return(compute_Xb(X, b))
    } else {
      return(compute_Xb(X^2, b))
    }
  }
}

compute_X2b_par = function(X, b, X_avg = 0) {
  if (is.list(X)) {
    n_var = unlist(lapply(X,get_ncol)) # number of variables for each element of X
    b_split = split_vector(b,n_var) # split b into a list of vectors
    return(
      foreach(i = 1:length(X), .combine = '+', .inorder = F) %dopar% {
        compute_X2b(X[[i]], b_split[[i]])
      }
    )
  } else {
    if (is.tfg_matrix(X)) {
      # X is boolean matrix, so X^2 = X
      return(compute_Xb(X, b))
    } else {
      return(compute_Xb(X^2, b))
    }
  }
}


# computes t((X - X_avg)^2) %*% y
compute_X2ty = function(X, y, X_avg = 0) {
  if (is.list(X)) {
    return(unlist(lapply(X, compute_X2ty, y = y)) - 2*compute_Xty(X, y)*X_avg + (X_avg^2 * sum(y)))
  } else {
    if (is.tfg_matrix(X)) {
      # X is boolean matrix, so X^2 = X
      return(as.numeric(compute_Xty(X, y)))
    } else {
      return(as.numeric(compute_Xty(X^2, y)))
    }
  }
}

compute_X2ty_par = function(X, y, X_avg = 0) {
  if (is.list(X)) {
    if (length(X_avg) != 1) {
      X_avg_split = split_vector(X_avg, sapply(X, get_ncol))
    } else {
      X_avg_split = lapply(1:length(X), function(i) X_avg)
    }
    return(
      foreach(i = 1:length(X), .combine = 'c') %dopar% {
        compute_X2ty(X[[i]], y) - (2 * compute_Xty(X[[i]], y) * X_avg_split[[i]])
        }  + (X_avg^2 * sum(y))
    )
  } else {
    if (is.tfg_matrix(X)) {
      # X is boolean matrix, so X^2 = X
      return(as.numeric(compute_Xty(X, y)))
    } else {
      return(as.numeric(compute_Xty(X^2, y)))
    }
  }
}


# weighted SER function, linear terms + stumps
# X is a list, first element corresponds to linear, others are stumps for variables
weighted_SER = function(X, Y, sigma2, init = list(V = NULL)) {
  if (length(sigma2) == 1) {
    sigma2 = rep(sigma2, length(Y))
  }
  
  inv_sigma2 = 1 / sigma2
  sum_inv_sigma2 = sum(inv_sigma2)
  w = inv_sigma2 / sum_inv_sigma2
  p = get_ncol(X)
  prior_weights = rep(1 / p, p)
  Y_avg = sum(Y * w)
  Y_cent = Y - Y_avg
  X_avg = compute_Xty(X, w) # vector of weighted avg of columns of X
  
  # tau_no_V = t(X_cent^2) %*% (1 / sigma2)
  #tau_no_V = compute_X2ty(X, inv_sigma2) - 2*compute_Xty(X, inv_sigma2)*X_avg + (X_avg^2 * sum_inv_sigma2)
  tau_no_V = compute_X2ty(X, inv_sigma2, X_avg)
  # nu = colSums(X_cent * Y_cent / sigma2) <=> t(X_cent) %*% (Y_cent / sigma2)
  nu = compute_Xty(X, Y_cent / sigma2) - (X_avg * sum(Y_cent / sigma2))
  
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
  intercept = as.numeric(Y_avg - sum(X_avg * beta_post_1))
  
  # mu1 = E[int + Xb] = E[Y_avg - X_avg'b + Xb]
  mu1 = intercept + compute_Xb(X, beta_post_1)
  # mu2 = E[(int + Xb)^2] = E[(Y_avg - X_avg'b + Xb)^2]
  #mu2 = Y_avg^2 + 2*Y_avg*(compute_Xb(X, beta_post_1) - sum(X_avg * beta_post_1)) + compute_X2b(X, beta_post_2) - 2*compute_Xb(X, beta_post_2*X_avg) + sum(X_avg^2 * beta_post_2)
  mu2 = Y_avg^2 + 2*Y_avg*(compute_Xb(X, beta_post_1) - sum(X_avg * beta_post_1)) + compute_X2b(X, beta_post_2, X_avg)
  
  KL_div = calc_KL(matrix(mu, ncol = 1), matrix(alpha, ncol = 1), matrix(sigma2_post, ncol = 1), V)
  
  return(list(mu1 = as.numeric(mu1), mu2 = as.numeric(mu2), KL_div = KL_div, alpha = alpha, mu = mu, sigma2_post = sigma2_post, intercept = intercept, V = V, X_avg = X_avg, Y_avg = Y_avg))
}

weighted_SER_par = function(X, Y, sigma2, init = list(V = NULL)) {
  if (length(sigma2) == 1) {
    sigma2 = rep(sigma2, length(Y))
  }
  
  inv_sigma2 = 1 / sigma2
  sum_inv_sigma2 = sum(inv_sigma2)
  w = inv_sigma2 / sum_inv_sigma2
  p = get_ncol(X)
  prior_weights = rep(1 / p, p)
  Y_avg = sum(Y * w)
  Y_cent = Y - Y_avg
  X_avg = compute_Xty_par(X, w) # vector of weighted avg of columns of X
  
  # tau_no_V = t(X_cent^2) %*% (1 / sigma2)
  #tau_no_V = compute_X2ty(X, inv_sigma2) - 2*compute_Xty(X, inv_sigma2)*X_avg + (X_avg^2 * sum_inv_sigma2)
  tau_no_V = compute_X2ty_par(X, inv_sigma2, X_avg)
  # nu = colSums(X_cent * Y_cent / sigma2) <=> t(X_cent) %*% (Y_cent / sigma2)
  nu = compute_Xty_par(X, Y_cent / sigma2) - (X_avg * sum(Y_cent / sigma2))
  
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
  intercept = as.numeric(Y_avg - sum(X_avg * beta_post_1))
  
  # mu1 = E[int + Xb] = E[Y_avg - X_avg'b + Xb]
  mu1 = intercept + compute_Xb_par(X, beta_post_1)
  # mu2 = E[(int + Xb)^2] = E[(Y_avg - X_avg'b + Xb)^2]
  #mu2 = Y_avg^2 + 2*Y_avg*(compute_Xb(X, beta_post_1) - sum(X_avg * beta_post_1)) + compute_X2b(X, beta_post_2) - 2*compute_Xb(X, beta_post_2*X_avg) + sum(X_avg^2 * beta_post_2)
  mu2 = Y_avg^2 + 2*Y_avg*(compute_Xb_par(X, beta_post_1) - sum(X_avg * beta_post_1)) +
    compute_X2b_par(X, beta_post_2, X_avg)
  
  KL_div = calc_KL(matrix(mu, ncol = 1), matrix(alpha, ncol = 1), matrix(sigma2_post, ncol = 1), V)
  
  return(list(mu1 = as.numeric(mu1), mu2 = as.numeric(mu2), KL_div = KL_div, alpha = alpha, mu = mu, sigma2_post = sigma2_post, intercept = intercept, V = V, X_avg = X_avg, Y_avg = Y_avg))
}

fitFn = function(X, Y, sigma2, init) {
  return(weighted_SER(X, Y, sigma2, init))
}

fitFnPar = function(X, Y, sigma2, init) {
  return(weighted_SER_par(X, Y, sigma2, init))
}

predFn = function(X_new, currentFit, moment = c(1, 2)) {
  beta_post_1 = currentFit$alpha * currentFit$mu
  if (moment == 1) {
    return(currentFit$intercept + compute_Xb(X_new, beta_post_1))
  } else if (moment == 2) {
    beta_post_2 = currentFit$alpha * (currentFit$sigma2_post + currentFit$mu^2)
    return(currentFit$Y_avg^2 + 2*currentFit$Y_avg*(compute_Xb(X_new, beta_post_1) - sum(currentFit$X_avg * beta_post_1)) + compute_X2b(X_new, beta_post_2, currentFit$X_avg))
  } else {
  } else {
    stop("`moment` must be either 1 or 2")
  }
}


# function to get weighted combination of fits
# returns (1 - rho)*fitOld + rho*fitNew
fitCombineFn = function(fitOld, fitNew, rho, X, Y, sigma2, predFn) {
  # find new X_avg and Y_avg (using all data points)
  inv_sigma2 = 1 / sigma2
  sum_inv_sigma2 = sum(inv_sigma2)
  w = inv_sigma2 / sum_inv_sigma2
  p = get_ncol(X)
  prior_weights = rep(1 / p, p)
  Y_avg = sum(Y * w)
  X_avg = compute_Xty(X, w)
  
  fitOut = list(
    alpha = (1 - rho)*fitOld$alpha + rho*fitNew$alpha, 
    mu = (1 - rho)*fitOld$mu + rho*fitNew$mu, 
    sigma2_post = (1 - rho)*fitOld$sigma2_post + rho*fitNew$sigma2_post, 
    intercept = (1 - rho)*fitOld$intercept + rho*fitNew$intercept, 
    V = (1 - rho)*fitOld$V + rho*fitNew$V, 
    X_avg = X_avg, Y_avg = Y_avg
  )
  
  fitOut$mu1 = predFn(X, fitOut, 1)
  fitOut$mu2 = predFn(X, fitOut, 2)
  fitOut$KL_div = calc_KL(matrix(fitOut$mu, ncol = 1), matrix(fitOut$alpha, ncol = 1), matrix(fitOut$sigma2_post, ncol = 1), fitOut$V)

  return(fitOut)
}


currentFitVector = function() {
  par = c(self$currentFit$alpha, self$currentFit$mu, self$currentFit$sigma2_post, self$currentFit$intercept, self$currentFit$V)
  return(par)
}

setCurrentFit = function(par) {
  p = get_ncol(self$X)
  par_remainder = par
  self$currentFit$alpha = par_remainder[1:p]
  par_remainder = par_remainder[-c(1:p)]
  self$currentFit$mu = par_remainder[1:p]
  par_remainder = par_remainder[-c(1:p)]
  self$currentFit$sigma2_post = par_remainder[1:p]
  par_remainder = par_remainder[-c(1:p)]
  self$currentFit$intercept = par_remainder[1]
  par_remainder = par_remainder[-1]
  self$currentFit$V = par_remainder[1]
  par_remainder = par_remainder[-1]
  return(par_remainder)
}
