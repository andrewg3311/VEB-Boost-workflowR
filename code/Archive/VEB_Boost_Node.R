### Node Object ###

VEBBoostNode <- R6::R6Class(
  "VEBBoostNode", 
  public = list(
    operator = NULL, # either "+" or "*" for internal nodes, NULL for terminal nodes
    
    fitFunction = NULL, # function that takes in predictors X, response Y, variances sigma2, and returns the fit
    # the fit must have fields mu1 (first moment), mu2 (second moment), KL_div (KL divergence from q to g),
    # and returns a function, predFunction(X_new, method = c(1, 2)) that takes in new data X and returns our prediction for the given moment
    
    currentFit = NULL, # current fit for fitting function
    
    ebCombineOperator = "+", # operator used to split current node. tries "+", then "*", then "." for locked
    
    nextNumber = 1,
    
    AddChildVEB = function(name, check = c("check", "no-warn", "no-check"), ...) { # add VEB node as child
      child = VEBBoostNode$new(as.character(name), check, ...)
      invisible(self$AddChildNode(child))
    },
    
    updateFit = function() { # function to update currentFit
      self$currentFit = self$fitFunction(X = self$X, Y = self$Y, sigma2 = self$sigma2, init = self$currentFit)
      invisible(self)
    }, 
    
    update_sigma2 = function() { # function to update sigma2
      self$sigma2 = ((sum(self$root$Y^2) - 2*sum(self$root$Y*self$root$mu1) + sum(self$root$mu2))) / length(self$root$Y)
      invisible(self)
    }, 
    
    convergeFit = function(tol = 1e-3, update_sigma2 = FALSE, update_ELBO_progress = FALSE, verbose = TRUE) {
      ELBOs = numeric(1000)
      ELBOs[1] = -Inf
      ELBOs[2] = self$root$ELBO
      i = 2
      while (abs(ELBOs[i] - ELBOs[i-1]) > tol) {
        self$root$Do(function(x) x$updateFit(), filterFun = function(x) x$isLeaf)
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
      
      self$fitFunction = NULL
      self$currentFit = NULL
      self$operator = operator
      self$name = combine_name
      
      self$AddChildNode(learner_copy)
      self$AddChildNode(learner)
      
      invisible(self$root) # return root, so we can assign entire tree to new value (in case we're splitting at root)
      
    }, 
    
    addLearner = function(learner_name, combine_name, fitFn = NULL, currentFit = NULL, tol = 1e-3) { # add a new learner to the tree
      base_learners = Traverse(self$root, filterFun = function(x) x$isLeaf & (x$KL_div > 0))
      var_inf = sapply(base_learners, function(x) x$variance_influence) # vector of var inf values
      #probs = var_inf / sum(var_inf)
      node_to_split = base_learners[[which.max(var_inf)]] # note: if all additive, this will lead to an unbalanced tree (so calculating moments, etc not as efficient)
      #node_to_split = base_learners[[sample.int(length(probs), size = 1, prob = probs)]] # randomly pick, weighted by influence
      #base_elbos = sapply(base_learners, function(x) x$ELBO)
      #node_to_split = base_learners[[which.min(base_elbos)]]
      
      node_clone = VEBBoostNode$new(node_to_split$name, fitFunction = node_to_split$fitFunction, currentFit = node_to_split$currentFit)
      node_clone$X = node_to_split$X
      node_clone$Y = node_to_split$Y
      node_clone$sigma2 = node_to_split$sigma2
      node_clone2 = node_clone$clone()
      
      if (is.null(fitFn)) {
        fitFn = node_clone$fitFunction
      }
      if (is.null(currentFit)) {
        # currentFitAdd = list(mu1 = rep(0, length(node_clone$Y)), mu2 = rep(0, length(node_clone$Y)), KL_div = 1e-10)
        # currentFitMult = list(mu1 = rep(1, length(node_clone$Y)), mu2 = rep(1, length(node_clone$Y)), KL_div = 1e-10)
        currentFitAdd = node_clone$currentFit
        currentFitMult = node_clone2$currentFit
      } else {
        currentFitAdd = currentFit
        currentFitMult = currentFit
      }
      
      # try adding nodes, see how we do
      add_node = VEBBoostNode$new(learner_name, fitFunction = fitFn, currentFit = currentFitAdd)
      add_learner = node_clone$AddSiblingVEB(add_node, "+", "add")
      add_learner$convergeFit(tol)
      
      # now try multiplying nodes, see how we do
      mult_node = VEBBoostNode$new(learner_name, fitFunction = fitFn, currentFit = currentFitMult)
      mult_learner = node_clone2$AddSiblingVEB(mult_node, "*", "mult")
      mult_learner$convergeFit(tol)
      
      if (add_learner$ELBO >= mult_learner$ELBO) {
        node_to_split$currentFit = add_learner[[node_to_split$name]]$currentFit
        new_node = VEBBoostNode$new(learner_name, fitFunction = add_node$fitFunction, currentFit = add_node$currentFit)
        node_to_split$AddSiblingVEB(new_node, "+", combine_name)
      } else {
        node_to_split$currentFit = mult_learner[[node_to_split$name]]$currentFit
        new_node = VEBBoostNode$new(learner_name, fitFunction = mult_node$fitFunction, currentFit = mult_node$currentFit)
        node_to_split$AddSiblingVEB(new_node, "*", combine_name)
      }
      
      invisible(self$root)
    },
    
    addLearnerAll = function(fitFn = NULL) { # to each leaf, add a "+" and "*"
      base_learners = Traverse(self$root, filterFun = function(x) x$isLeaf & (x$KL_div > 0))
      
      for (learner in base_learners) {
        fitFnLearner = fitFn
        if (is.null(fitFnLearner)) {
          fitFnLearner = learner$fitFunction
          name = learner$name
          add_learner_name = paste(name, "_add_learner", sep = "")
          mult_learner_name = paste(name, "_mult_learner", sep = "")
          add_combine_name = paste(name, "_add_combine", sep = "")
          mult_combine_name = paste(name, "_mult_combine", sep = "")
          add_fit = list(mu1 = rep(0, length(learner$Y)), mu2 = rep(0, length(learner$Y)), KL_div = 0)
          mult_fit = list(mu1 = rep(1, length(learner$Y)), mu2 = rep(1, length(learner$Y)), KL_div = 0)
          
          add_node = VEBBoostNode$new(add_learner_name, fitFunction = fitFnLearner, currentFit = add_fit)
          learner$AddSiblingVEB(add_node, "+", add_combine_name)
          
          mult_node = VEBBoostNode$new(mult_learner_name, fitFunction = fitFnLearner, currentFit = mult_fit)
          learner$children[[1]]$AddSiblingVEB(mult_node, "*", mult_combine_name)
        }
      }
      
      invisible(self$root)
    },
    
    addLearnerEB = function(learner_name, combine_name, fitFn = NULL, tol = 1e-3, V_tol = 1e-3) {
      if (self$ebCombineOperator == ".") { # if node is "locked"
        return(invisible(self))
      } else if (self$ebCombineOperator == "+") {
        eb_current_fit = list(mu1 = rep(0, length(self$Y)), mu2 = rep(0, length(self$Y)), KL_div = 0)
      } else {
        eb_current_fit = list(mu1 = rep(1, length(self$Y)), mu2 = rep(1, length(self$Y)), KL_div = 0)
      }
      
      node_clone = VEBBoostNode$new(self$name, fitFunction = self$fitFunction, currentFit = self$currentFit)
      node_clone$X = self$X
      node_clone$Y = self$Y
      node_clone$sigma2 = self$sigma2
      
      if (is.null(fitFn)) {
        fitFn = self$fitFunction
      }
      
      eb_node = VEBBoostNode$new(learner_name, fitFunction = fitFn, currentFit = eb_current_fit)
      eb_learner = node_clone$AddSiblingVEB(eb_node, self$ebCombineOperator, "eb_combine")
      eb_learner$convergeFit(tol)
      
      if (eb_learner$children[[2]]$currentFit$V < V_tol) { # if adding basically a constant
        if (self$ebCombineOperator == "+") { # try multiplicative combine
          self$ebCombineOperator = "*"
          return(self$addLearnerEB(learner_name, combine_name, fitFn, tol, V_tol))
        } else { # add multiplicative constant, "lock" node
          eb_node$fitFunction = fitFnConst
          eb_learner$convergeFit(tol)
          self$currentFit = eb_learner[[self$name]]$currentFit
          ebCombineOperator = self$ebCombineOperator
          self$ebCombineOperator = "."
          new_node = VEBBoostNode$new(learner_name, fitFunction = eb_node$fitFunction, currentFit = eb_node$currentFit, ebCombineOperator = ".")
          self$AddSiblingVEB(new_node, ebCombineOperator, combine_name)
        }
      } else {
        ebCombineOperator = self$ebCombineOperator
        self$ebCombineOperator = "+"
        self$currentFit = eb_learner[[self$name]]$currentFit
        new_node = VEBBoostNode$new(learner_name, fitFunction = eb_node$fitFunction, currentFit = eb_node$currentFit)
        self$AddSiblingVEB(new_node, ebCombineOperator, combine_name)
      }
      
      invisible(self)
    }, 
    
    convergeFitEB = function(tol = 1e-3, V_tol = 1e-3, update_sigma2 = FALSE, update_ELBO_progress = TRUE, verbose = TRUE) {
      base_learners = Traverse(self$root, filterFun = function(x) x$isLeaf & (x$ebCombineOperator != "."))
      for (learner in base_learners) { # change any near-constant leaf to a constant and seal it off
        if (learner$currentFit$V < V_tol) {
          learner$fitFunction = fitFnConst
          learner$ebCombineOperator = "."
          learner$updateFit()
        }
      }
      
      base_learners = Traverse(self$root, filterFun = function(x) x$isLeaf & (x$ebCombineOperator != "."))
      for (learner in base_learners) {
        learner_name = paste("mu_", learner$root$nextNumber, sep = '')
        combine_name = paste("combine_", learner$root$nextNumber, sep = '')
        learner$addLearnerEB(learner_name, combine_name)
        try({ # weird error, try to fix later, but still works
          learner$root$nextNumber = learner$root$nextNumber + 1
        }, silent = T)
        
        ELBOs = numeric(1000)
        ELBOs[1] = -Inf
        ELBOs[2] = self$root$ELBO
        i = 2
        while (abs(ELBOs[i] - ELBOs[i-1]) > tol) {
          self$root$Do(function(x) x$updateFit(), filterFun = function(x) x$isLeaf)
          if (update_sigma2) {
            self$root$update_sigma2()
          }
          i = i + 1
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
      }
      
      invisible(self$root)
    }, 
    
    predict.veb = function(X_new) { # function to get prediction on new data
      self$root$Do(function(node) {
        node$pred_mu1 = node$currentFit$predFunction(X_new, 1)
        node$pred_mu2 = node$currentFit$predFunction(X_new, 2)
        }, filterFun = function(x) x$isLeaf)
      invisible(self$root)
    }
    
  ),
  private = list(
    .X = NULL, # predictors for node (if NA, use value of parent)
    .Y = NA, # response, only not NA at root
    .sigma2 = NA, # variance, only not NA at root
    .ELBO_progress = list(-Inf), # list of ELBOs for each iteration of growing and fitting the tree
    .pred_mu1 = NULL, # prediction based on predFunction and given new data (first moment)
    .pred_mu2 = NULL # prediction based on predFunction and given new data (second moment)
  ), 
  
  active = list(

    mu1 = function(value) { # first moment of learner defined by sub-tree with this node as the root
      if (!missing(value)) {
        stop("`$mu1` cannot be modified directly", call. = FALSE)
      }
      
      if (self$isLeaf) {
        return(self$currentFit$mu1)
      }
      if (self$operator == "+") {
        return(apply(sapply(self$children, function(x) x$mu1), MARGIN = 1, sum))
      }
      if (self$operator == "*") {
        return(apply(sapply(self$children, function(x) x$mu1), MARGIN = 1, prod))
      }
    }, 
    
    mu2 = function(value) { # second moment of learner defined by sub-tree with this node as the root
      if (!missing(value)) {
        stop("`$mu2` cannot be modified directly", call. = FALSE)
      }
      
      if (self$isLeaf) {
        return(self$currentFit$mu2)
      }
      if (self$operator == "+") {
        return(apply(sapply(self$children, function(x) x$mu2), MARGIN = 1, sum) + 2*apply(sapply(self$children, function(x) x$mu1), MARGIN = 1, prod))
      }
      if (self$operator == "*") {
        return(apply(sapply(self$children, function(x) x$mu2), MARGIN = 1, prod))
      }
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
          return(private$.Y)
        }
        if (self$parent$operator == "+") {
          return(self$parent$Y - self$siblings[[1]]$mu1)
        }
        if (self$parent$operator == "*") {
          return(self$parent$Y * (self$siblings[[1]]$mu1 / self$siblings[[1]]$mu2))
        }
      } else {
        if (self$isRoot) {
          private$.Y = value
        } else {
          stop("`$Y` cannot be modified directly except at the root node", call. = FALSE)
        }
      }
    },
    
    sigma2 = function(value) { # variance for sub-tree
      if (missing(value)) {
        if (self$isRoot) {
          return(private$.sigma2)
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
      s2 = ifelse(length(self$sigma2) == 1, rep(self$sigma2, length(self$Y)), self$sigma2)
      return(
        (-.5 * sum(log(2*pi*s2))) - 
          (.5 * (sum(self$Y^2 / s2) - 2*sum(self$Y * self$mu1 / s2) + sum(self$mu2 / s2))) - 
          self$KL_div
      )
    }, 
    
    variance = function(value) { # variance of mean values at node
      if (!missing(value)) {
        stop("`$variance` cannot be modified directly", call. = FALSE)
      }
      
      return(self$mu2 - self$mu1^2)
    }, 
    
    # variance_partial_deriv = function(value) { # partial derivative of variance of mean w.r.t. variance at node
    #   if (!missing(value)) {
    #     stop("`$variance_partial_deriv` cannot be modified directly", call. = FALSE)
    #   }
    #   
    #   if (self$isRoot) {
    #     return(rep(1, length(self$Y)))
    #   }
    #   if (self$parent$operator == "+") {
    #     return(self$parent$variance_partial_deriv)
    #   } else {
    #     return(self$siblings[[1]]$variance * self$parent$variance_partial_deriv)
    #   }
    # }, 
    # 
    # variance_influence = function(value) { # weighted variance partial deriv, weighted by variance at node
    #   if (!missing(value)) {
    #     stop("`$variance_influence` cannot be modified directly", call. = FALSE)
    #   }
    #   
    #   if (all(self$Y == 0)) {
    #     return(-Inf)
    #   }
    #   return(weighted.mean(self$variance_partial_deriv, self$variance))
    #   # s2 = self$sigma2
    #   # if (length(self$sigma2) == 1) {
    #   #   s2 = rep(s2, length(self$Y))
    #   # }
    #   # return(mean(self$variance_partial_deriv * self$variance / s2))
    #   #return(sum(self$variance * self$variance_partial_deriv))
    # },
    
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
fitFnConst = function(X, Y, sigma2, init) {
  if (length(sigma2) == 1) {
    sigma2 = rep(sigma2, length(Y))
  }
  intercept = weighted.mean(Y, 1/sigma2)
  KL_div = 0
  predFunction = function(X_new, moment = c(1, 2)) {
    if (moment == 1) {
      return(rep(intercept, nrow(X_new)))
    } else if (moment == 2) {
      return(rep(intercept^2, nrow(X_new)))
    } else {
      stop("`moment` must be either 1 or 2")
    }
  }
  mu1 = predFunction(X, 1)
  mu2 = predFunction(X, 2)
  return(list(mu1 = mu1, mu2 = mu2, KL_div = KL_div, predFunction = predFunction))
}

# does take into account variability
# fitFnConst = function(X, Y, sigma2, init) {
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
# fitFnConst = function(X, Y, sigma2, init) {
#   if (length(sigma2) == 1) {
#     sigma2 = rep(sigma2, length(Y))
#   }
#   intercept = weighted.mean(Y, 1/sigma2)
#   KL_div = .5 * (log(sum(1 / sigma2)) - 1)
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
