### Node Object ###

### This version is more computationally efficient by storing the moments in the internal nodes as well, rather than having them be actives.
### When we update a base learner's fit, we pass the updates up the tree and update the moments in the internal nodes.
### This should be faster computationally, but more memory intensive, since we have to store more data.

### THIS version modifies updateMoments() and updateFit() so that you only update what's necessary, hopefully decreasing over-head

### THIS version stores $mu1 and $mu2 as big.matrix objects, so that we can update nodes in parallel w/out having to copy everything to each node
### Each node's private$.mu1 and private$.mu2 stores a description for the node's mu1 and mu2 big.matrix
### The active methods $mu1 and $mu2 first attach the matrix, then either return or modify the values therein

VEBBoostNodeBM <- R6::R6Class(
  "VEBBoostNodeBM", 
  public = list(
    operator = NULL, # either "+" or "*" for internal nodes, NULL for terminal nodes
    
    fitFunction = NULL, # function that takes in predictors X, response Y, variances sigma2, and returns the fit
    # the fit must have fields mu1 (first moment), mu2 (second moment), KL_div (KL divergence from q to g),
    # and returns a function, predFunction(X_new, method = c(1, 2)) that takes in new data X and returns our prediction for the given moment
    
    currentFit = NULL, # current fit for fitting function
    
    predFunction = NULL, # function to predict based on current fit
    
    family = "gaussian",
    
    AddChildVEB = function(name, check = c("check", "no-warn", "no-check"), ...) { # add VEB node as child
      child = VEBBoostNodeBM$new(as.character(name), check, ...)
      invisible(self$AddChildNode(child))
    },
    
    updateMoments = function() { # after updating the fit, pass changes to moments up to parent node
      if (!self$isLeaf) { # if not at a leaf, update moments
        if (self$operator == "+") {
          children_mu1 = sapply(self$children, function(x) x$mu1)
          self$mu1 <- as.numeric(apply(children_mu1, MARGIN = 1, sum))
          self$mu2 <- as.numeric(apply(sapply(self$children, function(x) x$mu2), MARGIN = 1, sum) + 2*apply(children_mu1, MARGIN = 1, prod))
        } else {
          self$mu1 <- as.numeric(apply(sapply(self$children, function(x) x$mu1), MARGIN = 1, prod))
          self$mu2 <- as.numeric(apply(sapply(self$children, function(x) x$mu2), MARGIN = 1, prod))
        }
      }
      
      return(invisible(self))
    }, 
    
    updateMomentsAll = function() { # after updating the fit, pass changes to moments up to internal nodes
      if (!self$isLeaf) { # if not at a leaf, update moments
        if (self$operator == "+") {
          children_mu1 = sapply(self$children, function(x) x$mu1)
          self$mu1 <- as.numeric(apply(children_mu1, MARGIN = 1, sum))
          self$mu2 <- as.numeric(apply(sapply(self$children, function(x) x$mu2), MARGIN = 1, sum) + 2*apply(children_mu1, MARGIN = 1, prod))
        } else {
          self$mu1 <- as.numeric(apply(sapply(self$children, function(x) x$mu1), MARGIN = 1, prod))
          self$mu2 <- as.numeric(apply(sapply(self$children, function(x) x$mu2), MARGIN = 1, prod))
        }
      }
      if (self$isRoot) { # if at root, stop
        return(invisible(self))
      } else { # else, update parents moments
        return(invisible(self$parent$updateMomentsAll()))
      }
    }, 
    
    updateFit = function() { # function to update currentFit
      if (self$isLeaf) {
        self$currentFit = self$fitFunction(X = self$X, Y = self$Y, sigma2 = self$sigma2, init = self$currentFit)
        self$mu1 <- self$predFunction(self$X, self$currentFit, 1)
        self$mu2 <- self$predFunction(self$X, self$currentFit, 2)
        self$updateMomentsAll()
      }
      invisible(self)
    }, 
    
    update_sigma2 = function() { # function to update sigma2
      self$sigma2 = ((sum(self$root$Y^2) - 2*sum(self$root$Y*self$root$mu1) + sum(self$root$mu2))) / length(self$root$Y)
      invisible(self)
    }, 
    
    convergeFit = function(tol = 1e-3, update_sigma2 = FALSE, update_ELBO_progress = TRUE, verbose = TRUE) {
      ELBOs = numeric(1000)
      ELBOs[1] = -Inf
      ELBOs[2] = self$root$ELBO
      i = 2
      base_learners = Traverse(self$root, filterFun = function(x) x$isLeaf, traversal = 'post-order')
      while (abs(ELBOs[i] - ELBOs[i-1]) > tol) {
        #self$root$Do(function(x) try({x$updateFit()}, silent = T), traversal = 'post-order')
        newFits <- foreach(learner = base_learners, .packages = c('bigmemory', 'data.tree')) %dopar% {
          learner$updateFit()
          learner$currentFit
        }
        
        for (j in 1:length(base_learners)) {
          base_learners[[j]]$currentFit <- newFits[[j]]
        }
        # need the try to avoid error/bug in data.tree or R6, that random bug where is says trying to apply non-function
        if (update_sigma2) {
          self$root$update_sigma2()
        }
        i = i+1
        if (i > length(ELBOs)) { # double size of ELBOs for efficiency rather than making it bigger each iteration
          ELBOs = c(ELBOs, rep(0, i))
        }
        ELBOs[i] = self$root$ELBO
        if (verbose & ((i %% 100) == 0)) {
          cat(paste("ELBO: ", ELBOs[i], sep = ""))
          cat("\n")
        }
      }
      ELBOs = ELBOs[2:i]
      
      if (update_ELBO_progress) {
        self$ELBO_progress = ELBOs
      }
      invisible(self$root)
    }, 
    
    AddSiblingVEB = function(learner, operator = c("+", "*"), combine_name) { # function to add subtree as sibling to given node, combining with operator
      # learner is tree to add to self
      # operator is how to combine them
      # name is what to call the combining node
      
      learner_copy = self$clone()
      learner_copy$X = NULL
      #learner_copy$mu = bigmemory::describe(bigmemory::deepcopy(self$mu))
      #eval(parse(text = paste(combine_name, "_mu", " <- ", sep = "")))
      
      # make copy in storage environment
      assign(paste(combine_name, "_mu", sep = ""), bigmemory::deepcopy(self$mu), envir = eval(parse(text = self$storageEnvironment)))
      # assign pointer to copy
      eval(parse(text = paste("learner_copy$mu <- bigmemory::describe(", self$storageEnvironment, "$", combine_name, "_mu)", sep = "")))
      
      self$fitFunction = NULL
      self$currentFit = NULL
      self$predFunction = NULL
      self$operator = operator
      self$name = combine_name
      
      self$AddChildNode(learner_copy)
      self$AddChildNode(learner)
      
      learner_copy$parent$updateMoments()
      
      invisible(self$root) # return root, so we can assign entire tree to new value (in case we're splitting at root)
      
    }, 
    
    lockLearners = function(V_tol = 1e-3) { # lock learners that should be locked
      base_learners = Traverse(self$root, traversal = 'post-order', filterFun = function(x) x$isLeaf & !x$isLocked)
      for (learner in base_learners) { # change any near-constant leaf to a constant and seal it off
        if (learner$currentFit$V < V_tol) {
          learner$isLocked = TRUE
          learner$fitFunction = fitFnConstComp
          learner$predFunction = predFnConstComp
          learner$updateFit()
          #try({learner$updateMomentsAll()}, silent = T) # needed to avoid "attempt to apply non-function" error
        }
      }
      
      base_learners = Traverse(self$root, filterFun = function(x) x$isLeaf & !x$isLocked)
      for (learner in base_learners) {
        if (learner$isRoot) { # if root, not locked, do this to avoid errors in next if statement
          next
        }
        if ((learner$siblings[[1]]$isLocked) && learner$parent$siblings[[1]]$isLocked) {
          learner$isLocked = TRUE
        }
      }
      
      return(invisible(self$root))
    }, 
    
    addLearnerAll = function(V_tol = 1e-3) { # to each leaf, add a "+" and "*"
      self$root$lockLearners(V_tol) # lock learners
      
      base_learners = Traverse(self$root, filterFun = function(x) x$isLeaf & !x$isLocked)
      for (learner in base_learners) {
        fitFn = learner$fitFunction
        predFn = learner$predFunction
        learner_name = paste("mu_", learner$root$leafCount, sep = '')
        combine_name = paste("combine_", learner$root$leafCount, sep = '')
        
        add_node = VEBBoostNodeBM$new(learner_name, fitFunction = fitFn, predFunction = predFn, currentFit = list(KL_div = 0))
        assign(paste(learner_name, "_mu", sep = ""), bigmemory::deepcopy(learner$mu), envir = eval(parse(text = learner$storageEnvironment)))
        # assign pointer to copy
        eval(parse(text = paste("add_node$mu <- bigmemory::describe(", learner$storageEnvironment, "$", learner_name, "_mu)", sep = "")))
        #add_node$mu <- bigmemory::describe(bigmemory::deepcopy(learner$mu))
        add_node$mu <- 0
        learner$AddSiblingVEB(add_node, "+", combine_name)
        
        learner_name = paste("mu_", learner$root$leafCount, sep = '')
        combine_name = paste("combine_", learner$root$leafCount, sep = '')
        
        mult_node = VEBBoostNodeBM$new(learner_name, fitFunction = fitFn, predFunction = predFn, currentFit = list(KL_div = 0))
        assign(paste(learner_name, "_mu", sep = ""), bigmemory::deepcopy(learner$mu), envir = eval(parse(text = learner$storageEnvironment)))
        # assign pointer to copy
        eval(parse(text = paste("mult_node$mu <- bigmemory::describe(", learner$storageEnvironment, "$", learner_name, "_mu)", sep = "")))
        #mult_node$mu <- bigmemory::describe(bigmemory::deepcopy(learner$mu))
        mult_node$mu <- 1
        learner$children[[1]]$AddSiblingVEB(mult_node, "*", combine_name)
      }
      
      invisible(self$root)
    },
    
    convergeFitAll = function(tol = 1e-3, V_tol = 1e-3, update_sigma2 = FALSE, update_ELBO_progress = TRUE, verbose = TRUE) {
      self$addLearnerAll(V_tol)
      #self$Do(function(node) node$updateMoments(), traversal = 'post-order')
      self$convergeFit(tol, update_sigma2, update_ELBO_progress, verbose)
      return(invisible(self$root))
    },
    
    predict.veb = function(X_new, moment = c(1, 2)) { # function to get prediction on new data
      self$root$Do(function(node) {
        if (1 %in% moment) {
          node$pred_mu1 = node$predFunction(X_new, node$currentFit, 1)
        }
        if (2 %in% moment) {
          node$pred_mu2 = node$predFunction(X_new, node$currentFit, 2)
        }
      }, filterFun = function(x) x$isLeaf)
      invisible(self$root)
    }
    
  ),
  private = list(
    .X = NULL, # predictors for node (if NA, use value of parent)
    .Y = NA, # response, only not NA at root
    .sigma2 = NA, # variance, only not NA at root
    .mu = NA, # description for matrix w/ 1st column current first moments, 2nd column current second moments
    .storageEnvironment = NULL, # string of environment for where all mu1 and mu2 big.matrices are stored
    .ELBO_progress = list(-Inf), # list of ELBOs for each iteration of growing and fitting the tree
    .pred_mu1 = NULL, # prediction based on predFunction and given new data (first moment)
    .pred_mu2 = NULL, # prediction based on predFunction and given new data (second moment)
    .isLocked = FALSE, # locked <=> V < V_tol, or both learners directly connected (sibling and parent's sibling) are constant
    #.ensemble = NULL, # only really used in mult-class case, otherwise will refer to root
    .alpha = 0 # used in multi-class learner
  ), 
  
  active = list(
    
    ensemble = function(value) {
      if (missing(value)) {
        if (self$isRoot) {
          if (identical(parent.env(self), emptyenv())) {
            return(invisible(self))
          } else {
            return(invisible(parent.env(self)))
          }
        } else {
          return(invisible(self$root$ensemble))
        }
      }
      
      if (!self$isRoot) {
        stop("`$ensemble' cannot be modified except at the root", call. = FALSE)
      } else {
        parent.env(self) = value
      }
    },

    isEnsemble = function(value) { # TRUE if is ensemble (root AND not part of higher ensemble), else FALSE
      if (!missing(value)) {
        stop("`$isEnsemble' cannot be modified directly", call. = FALSE)
      }
      if (self$isRoot) {
        if (identical(parent.env(self), emptyenv())) {
          return(TRUE)
        }
      }
      return(FALSE)
    },
    
    storageEnvironment = function(value) {
      if (!missing(value)) {
        if (is.environment(eval(parse(text = value)))) {
          private$.storageEnvironment = value
        } else {
          stop(paste("No environment named ", value, " found", sep = ""))
        }
      } else {
        if (!is.null(private$.storageEnvironment)) {
          return(private$.storageEnvironment)
        } else {
          if (self$isRoot) {
            return(".GlobalEnv")
          } else {
            return(self$parent$storageEnvironment)
          }
        }
      }
    }, 
    
    mu = function(value) {
      if (!missing(value)) {
        if (class(value) == "big.matrix.descriptor") {
          private$.mu <- value
        } else {
          mu_bm <- bigmemory::attach.big.matrix(private$.mu)
          mu_bm[] <- value
        }
        return(invisible(self))
      } else {
        mu_bm <- bigmemory::attach.big.matrix(private$.mu)
        return(mu_bm[])
      }
    }, 
    
    mu1 = function(value) {
      mu_bm <- bigmemory::attach.big.matrix(private$.mu)
      if (!missing(value)) {
        mu_bm[, 1] <- value
        return(invisible(self))
      } else {
        return(mu_bm[, 1])
      }
    }, 
    
    mu2 = function(value) {
      mu_bm <- bigmemory::attach.big.matrix(private$.mu)
      if (!missing(value)) {
        mu_bm[, 2] <- value
        return(invisible(self))
      } else {
        return(mu_bm[, 2])
      }
    }, 
    
    mu1_div_mu2 = function(value) {
      if (!missing(value)) {
        stop("`$mu1_div_mu21 cannot be modified directly")
      }
      mu_bm <- bigmemory::attach.big.matrix(private$.mu)
      return(mu_bm[, 1] / mu_bm[, 2])
    }, 
    
    KL_div = function(value) { # KL divergence from q to g of learner defined by sub-tree with this node as the root
      if (!missing(value)) {
        stop("`$KL_div` cannot be modified directly", call. = FALSE)
      }
      
      if (self$isLeaf) {
        return(self$currentFit$KL_div)
      }
      return(sum(sapply(self$children, function(x) x$KL_div)))
    }, 
    
    X = function(value) { # predictor for node
      if (missing(value)) {
        if (self$isRoot) {
          if (self$isEnsemble) {
            return(private$.X)
          } else {
            return(self$ensemble$X)
          }
        }
        if (is.null(private$.X)) {
          return(self$parent$X)
        } else {
          return(private$.X)
        }
      }
      private$.X = value
    }, 
    
    Y = function(value) { # response for sub-tree
      if (missing(value)) {
        if (self$isRoot) {
          if (self$root$family == "gaussian") {
            return(private$.Y)
          } else if (self$root$family == "binomial") {
            d = self$root$d
            return((private$.Y - .5 + (self$alpha * d)) / d) # alpha*d needed for multi-class case
          } else {
            stop("family must be either 'gaussian' or 'binomial")
          }
        }
        if (self$parent$operator == "+") {
          return(self$parent$Y - self$siblings[[1]]$mu1)
        }
        if (self$parent$operator == "*") {
          return(self$parent$Y * self$siblings[[1]]$mu1_div_mu2)
        }
      } else {
        if (self$isRoot) {
          private$.Y = value
        } else {
          stop("`$Y` cannot be modified directly except at the root node", call. = FALSE)
        }
      }
    },
    
    raw_Y = function(value) { # raw value of private$.Y, only needed for logistic ELBO
      if (!missing(value)) {
        stop("`$raw_Y` cannot be modified directly", call. = FALSE)
      }
      if (self$isRoot) {
        return(private$.Y)
      } else {
        return(self$parent$raw_Y)
      }
    }, 
    
    sigma2 = function(value) { # variance for sub-tree
      if (missing(value)) {
        if (self$isRoot) {
          if (self$root$family == "gaussian") {
            return(private$.sigma2)
          } else if (self$root$family == "binomial") {
            return(1 / self$root$d)
          } else {
            stop("family must be either 'gaussian' or 'binomial")
          }
        }
        if (self$parent$operator == "+") {
          return(self$parent$sigma2)
        }
        if (self$parent$operator == "*") {
          return(self$parent$sigma2 / self$siblings[[1]]$mu2)
        }
      } else {
        if (self$isRoot) {
          private$.sigma2 = value
        } else {
          stop("`$sigma2` cannot be modified directly except at the root node", call. = FALSE)
        }
      } 
    },
    
    ELBO_progress = function(value) { # ELBO progress of tree
      if (missing(value)) {
        if (self$isRoot) {
          return(private$.ELBO_progress)
        } else {
          return(self$root$ELBO_progress)
        }
      } else {
        if (self$isRoot) { # append ELBOs in list
          private$.ELBO_progress[[length(private$.ELBO_progress) + 1]] = value
        } else {
          stop("`$ELBO_progress` cannot be modified directly except at the root node", call. = FALSE)
        }
      }
    }, 
    
    ELBO = function(value) { # ELBO for entire tree
      if (!missing(value)) {
        stop("`$ELBO` cannot be modified directly", call. = FALSE)
      }
      # make s2 vector of variances
      #s2 = ifelse(length(self$sigma2) == 1, rep(self$sigma2, length(self$Y)), self$sigma2) # something weird w/ this if-else, not sure why
      if (self$root$family == "gaussian") {
        s2 = self$sigma2
        if (length(s2) == 1) {
          s2 = rep(s2, length(self$Y))
        }
        return(
          (-.5 * sum(log(2*pi*s2))) - 
            (.5 * (sum(self$Y^2 / s2) - 2*sum(self$Y * self$mu1 / s2) + sum(self$mu2 / s2))) - 
            self$KL_div
        )
      } else if (self$root$family == "binomial") {
        d = self$root$d
        xi = self$root$xi
        return(
          sum(log(g(xi))) + sum((xi / 2)*(d*xi - 1)) + sum((self$root$raw_Y - .5) * self$root$mu1) - 
            .5*sum(self$root$mu2 * d) - self$KL_div
        )
      } else {
        stop("family must be either 'gaussian' or 'binomial")
      }
    }, 
    
    alpha = function(value) { # used in multi-class learner
      if (!missing(value)) {
        stop("`$alpha' cannot be modified directly except at the ensemble level", call. = FALSE)
      }
      if (self$isEnsemble) { # if isEnsemble, e.g. if we're in linear or logistic case, return root$.alpha (SHOULD BE 0 IN THIS CASE)
        return(private$.alpha)
      }
      return(self$ensemble$alpha) # otherwise, return ensemble's alpha
    }, 
    
    xi = function(value) { # optimal variational parameters, set to +sqrt(mu2)
      if (!missing(value)) {
        stop("`$xi` cannot be modified directly", call. = FALSE)
      }
      return(sqrt(self$root$mu2 + self$alpha^2 - 2*self$alpha*self$root$mu1))
    }, 
    
    d = function(value) { # d == 1/xi * (g(xi) - .5), n x K matrix
      if (!missing(value)) {
        stop("'$d' cannot be modified directly", call. = FALSE)
      }
      xi = self$xi
      g_xi = g(xi) # matrix of g(xi_i,k), pre-compute once
      d = ((g_xi - .5) / xi) # matrix of (g(xi_i,k) - .5) / xi_i,k, pre-compute once
      d[xi == 0] = .25 # case of 0/0 (e.g. x_i is all 0), use L'Hopital
      return(d)
    }, 
    
    isLocked = function(value) { # locked <=> V < V_tol, or both learners directly connected (sibling and parent's sibling) are constant
      if (missing(value)) {
        if (self$isLeaf) {
          return(private$.isLocked)
        } else {
          return(all(sapply(self$children, function(x) x$isLocked)))
        }
      } else {
        if (self$isLeaf) {
          private$.isLocked = value
        } else {
          stop("`$isLocked` cannot be modified directly except at leaf nodes", call. = FALSE)
        }
      }
    },
    
    
    variance = function(value) { # variance of mean values at node
      if (!missing(value)) {
        stop("`$variance` cannot be modified directly", call. = FALSE)
      }
      
      return(self$mu2 - self$mu1^2)
    }, 
    
    pred_mu1 = function(value) { # predicted first moment given new data
      if (!missing(value)) {
        if (self$isLeaf) {
          private$.pred_mu1 = value
          return(invisible(self))
        } else {
          stop("`$pred_mu1` cannot be modified directly except at leaf nodes by calling `$predict.veb`", call. = FALSE)
        }
      }
      if (self$isLeaf) {
        return(private$.pred_mu1)
      } else if (self$operator == "+") {
        return(apply(sapply(self$children, function(x) x$pred_mu1), MARGIN = 1, sum))
      } else {
        return(apply(sapply(self$children, function(x) x$pred_mu1), MARGIN = 1, prod))
      }
    },
    
    pred_mu2 = function(value) { # predicted second moment given new data
      if (!missing(value)) {
        if (self$isLeaf) {
          private$.pred_mu2 = value
          return(invisible(self))
        } else {
          stop("`$pred_mu2` cannot be modified directly except at leaf nodes by calling `$predict.veb`", call. = FALSE)
        }
      }
      if (self$isLeaf) {
        return(private$.pred_mu2)
      } else if (self$operator == "+") {
        return(apply(sapply(self$children, function(x) x$pred_mu2), MARGIN = 1, sum) + 2*apply(sapply(self$children, function(x) x$pred_mu1), MARGIN = 1, prod))
      } else {
        return(apply(sapply(self$children, function(x) x$pred_mu2), MARGIN = 1, prod))
      }
    }
  ), 
  inherit = data.tree::Node, 
  lock_objects = FALSE
)

# does not take into account variability
fitFnConstComp = function(X, Y, sigma2, init) {
  if (length(sigma2) == 1) {
    sigma2 = rep(sigma2, length(Y))
  }
  intercept = weighted.mean(Y, 1/sigma2)
  KL_div = 0
  
  return(list(intercept = intercept, KL_div = KL_div))
}

predFnConstComp = function(X_new, currentFit, moment = c(1, 2)) {
  if (moment == 1) {
    return(rep(currentFit$intercept, get_nrow(X_new)))
  } else if (moment == 2) {
    return(rep(currentFit$intercept^2, get_nrow(X_new)))
  } else {
    stop("`moment` must be either 1 or 2")
  }
}



# does take into account variability
# fitFnConstComp = function(X, Y, sigma2, init) {
#   if (length(sigma2) == 1) {
#     sigma2 = rep(sigma2, length(Y))
#   }
#   intercept = weighted.mean(Y, 1/sigma2)
#   KL_div = 0
#   predFunction = function(X_new, moment = c(1, 2)) {
#     if (moment == 1) {
#       return(rep(intercept, nrow(X_new)))
#     } else if (moment == 2) {
#       return(rep(intercept^2 + sum(1 / sigma2)^(-1), nrow(X_new)))
#     } else {
#       stop("`moment` must be either 1 or 2")
#     }
#   }
#   mu1 = predFunction(X, 1)
#   mu2 = predFunction(X, 2)
#   return(list(mu1 = mu1, mu2 = mu2, KL_div = KL_div, predFunction = predFunction))
# }

# more appropriate KL divergence (still excludes constant)
# fitFnConstComp = function(X, Y, sigma2, init) {
#   if (length(sigma2) == 1) {
#     sigma2 = rep(sigma2, length(Y))
#   }
#   intercept = weighted.mean(Y, 1/sigma2)
#   KL_div = .5 * (log(sum(1 / sigma2)) - 1)
#   mu1 = rep(intercept, nrow(X))
#   mu2 = rep(intercept^2 + sum(1 / sigma2)^(-1), nrow(X))
#   return(list(mu1 = mu1, mu2 = mu2, KL_div = KL_div))
# }
# 
# predFnConstComp = function(X_new, currentFit, moment = c(1, 2)) {
#   if (moment == 1) {
#     return(rep(currentFit$intercept, nrow(X_new)))
#   } else if (moment == 2) {
#     return(rep(currentFit$intercept^2 + sum(1 / sigma2)^(-1), nrow(X_new)))
#   } else {
#     stop("`moment` must be either 1 or 2")
#   }
# }

# does EB for variance of intercept
# fitFnConstComp = function(X, Y, sigma2, init) {
#   if (length(sigma2) == 1) {
#     sigma2 = rep(sigma2, length(Y))
#   }
#   sum_inv_sigma2 = sum(1 / sigma2)
#   nu_int = sum(Y / sigma2)
#   V_int = max(0, ((nu_int^2 / sum_inv_sigma2) - 1) / sum_inv_sigma2)
#   tau_int = (1 / V_int) + sum_inv_sigma2
#   mu1 = ifelse(V_int == 0, 0, nu_int / tau_int)
#   mu2 = ifelse(V_int == 0, 0, 1 / tau_int + mu1^2)
#   KL_div = ifelse(V_int == 0, 0, ((log(V_int) + log(tau_int)) + (((1 / tau_int) + (mu1^2)) /  V_int) - 1) / 2)
#   return(list(mu1 = rep(mu1, length(Y)), mu2 = rep(mu2,length(Y)), KL_div = KL_div, V_int = V_int))
# }
# 
# predFnConstComp = function(X_new, currentFit, moment = c(1, 2)) {
#   if (moment == 1) {
#     return(currentFit$mu1)
#   } else if (moment == 2) {
#     return(currentFit$mu2)
#   } else {
#     stop("`moment` must be either 1 or 2")
#   }
# }

g = function(x) {
  1 / (1 + exp(-x))
}

