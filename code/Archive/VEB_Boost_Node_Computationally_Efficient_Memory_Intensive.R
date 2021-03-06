### Node Object ###

### This version is more computationally efficient by storing the moments in the internal nodes as well, rather than having them be actives.
### When we update a base learner's fit, we pass the updates up the tree and update the moments in the internal nodes.
### This should be faster computationally, but more memory intensive, since we have to store more data.

VEBBoostNodeComp <- R6::R6Class(
  "VEBBoostNodeComp", 
  public = list(
    operator = NULL, # either "+" or "*" for internal nodes, NULL for terminal nodes
    
    fitFunction = NULL, # function that takes in predictors X, response Y, variances sigma2, and returns the fit
    # the fit must have fields mu1 (first moment), mu2 (second moment), KL_div (KL divergence from q to g),
    # and returns a function, predFunction(X_new, method = c(1, 2)) that takes in new data X and returns our prediction for the given moment
    
    currentFit = NULL, # current fit for fitting function
    
    predFunction = NULL, # function to predict based on current fit
    
    family = "gaussian",
    
    ebCombineOperator = "+", # operator used to split current node. tries "+", then "*", then "." for locked
    
    nextNumber = 1,
    
    AddChildVEB = function(name, check = c("check", "no-warn", "no-check"), ...) { # add VEB node as child
      child = VEBBoostNodeComp$new(as.character(name), check, ...)
      invisible(self$AddChildNode(child))
    },
    
    updateMoments = function() { # after updating the fit, pass changes to moments up to internal nodes
      if (!self$isLeaf) { # if not at a leaf, update moments
        if (self$operator == "+") {
          private$.mu1 = as.numeric(apply(sapply(self$children, function(x) x$mu1), MARGIN = 1, sum))
          private$.mu2 = as.numeric(apply(sapply(self$children, function(x) x$mu2), MARGIN = 1, sum) + 2*apply(sapply(self$children, function(x) x$mu1), MARGIN = 1, prod))
        } else {
          private$.mu1 = as.numeric(apply(sapply(self$children, function(x) x$mu1), MARGIN = 1, prod))
          private$.mu2 = as.numeric(apply(sapply(self$children, function(x) x$mu2), MARGIN = 1, prod))
        }
      }
      if (self$isRoot) { # if at root, stop
        return(invisible(self))
      } else { # else, update parents moments
        return(invisible(self$parent$updateMoments()))
      }
    }, 
    
    updateFit = function() { # function to update currentFit
      self$currentFit = self$fitFunction(X = self$X, Y = self$Y, sigma2 = self$sigma2, init = self$currentFit)
      self$updateMoments()
      invisible(self)
    }, 
    
    update_sigma2 = function() { # function to update sigma2
      self$sigma2 = ((sum(self$root$Y^2) - 2*sum(self$root$Y*self$root$mu1) + sum(self$root$mu2))) / length(self$root$Y)
      invisible(self)
    }, 
    
    # update_sigma2 = function(stepSize = .01) { # function to update sigma2
    #   # stepSize is how much weight to give new estimate
    #   sigma2_new = ((sum(self$root$Y^2) - 2*sum(self$root$Y*self$root$mu1) + sum(self$root$mu2))) / length(self$root$Y)
    #   self$sigma2 = (stepSize * sigma2_new) + ((1 - stepSize) * self$sigma2)
    #   invisible(self)
    # }, 
    
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
      self$predFunction = NULL
      self$operator = operator
      self$name = combine_name
      
      self$AddChildNode(learner_copy)
      self$AddChildNode(learner)
      
      learner_copy$updateMoments()
      
      invisible(self$root) # return root, so we can assign entire tree to new value (in case we're splitting at root)
      
    }, 
    
    # addLearner = function(learner_name, combine_name, fitFn = NULL, currentFit = NULL, predFn, tol = 1e-3) { # add a new learner to the tree
    #   base_learners = Traverse(self$root, filterFun = function(x) x$isLeaf & (x$KL_div > 0))
    #   var_inf = sapply(base_learners, function(x) x$variance_influence) # vector of var inf values
    #   #probs = var_inf / sum(var_inf)
    #   node_to_split = base_learners[[which.max(var_inf)]] # note: if all additive, this will lead to an unbalanced tree (so calculating moments, etc not as efficient)
    #   #node_to_split = base_learners[[sample.int(length(probs), size = 1, prob = probs)]] # randomly pick, weighted by influence
    #   #base_elbos = sapply(base_learners, function(x) x$ELBO)
    #   #node_to_split = base_learners[[which.min(base_elbos)]]
    #   
    #   node_clone = VEBBoostNodeComp$new(node_to_split$name, fitFunction = node_to_split$fitFunction, currentFit = node_to_split$currentFit)
    #   node_clone$X = node_to_split$X
    #   node_clone$Y = node_to_split$Y
    #   node_clone$sigma2 = node_to_split$sigma2
    #   node_clone2 = node_clone$clone()
    #   
    #   if (is.null(fitFn)) {
    #     fitFn = node_clone$fitFunction
    #   }
    #   if (is.null(predFn)) {
    #     predFn = node_clone$predFunction
    #   }
    #   if (is.null(currentFit)) {
    #     # currentFitAdd = list(mu1 = rep(0, length(node_clone$Y)), mu2 = rep(0, length(node_clone$Y)), KL_div = 1e-10)
    #     # currentFitMult = list(mu1 = rep(1, length(node_clone$Y)), mu2 = rep(1, length(node_clone$Y)), KL_div = 1e-10)
    #     currentFitAdd = node_clone$currentFit
    #     currentFitMult = node_clone2$currentFit
    #   } else {
    #     currentFitAdd = currentFit
    #     currentFitMult = currentFit
    #   }
    #   
    #   # try adding nodes, see how we do
    #   add_node = VEBBoostNodeComp$new(learner_name, fitFunction = fitFn, currentFit = currentFitAdd)
    #   add_learner = node_clone$AddSiblingVEB(add_node, "+", "add")
    #   add_learner$convergeFit(tol)
    #   
    #   # now try multiplying nodes, see how we do
    #   mult_node = VEBBoostNodeComp$new(learner_name, fitFunction = fitFn, currentFit = currentFitMult)
    #   mult_learner = node_clone2$AddSiblingVEB(mult_node, "*", "mult")
    #   mult_learner$convergeFit(tol)
    #   
    #   if (add_learner$ELBO >= mult_learner$ELBO) {
    #     node_to_split$currentFit = add_learner[[node_to_split$name]]$currentFit
    #     new_node = VEBBoostNodeComp$new(learner_name, fitFunction = add_node$fitFunction, currentFit = add_node$currentFit)
    #     node_to_split$AddSiblingVEB(new_node, "+", combine_name)
    #   } else {
    #     node_to_split$currentFit = mult_learner[[node_to_split$name]]$currentFit
    #     new_node = VEBBoostNodeComp$new(learner_name, fitFunction = mult_node$fitFunction, currentFit = mult_node$currentFit)
    #     node_to_split$AddSiblingVEB(new_node, "*", combine_name)
    #   }
    #   
    #   invisible(self$root)
    # },
    
    lockLearners = function(V_tol = 1e-3) { # lock learners that should be locked
      base_learners = Traverse(self$root, filterFun = function(x) x$isLeaf & !x$isLocked)
      for (learner in base_learners) { # change any near-constant leaf to a constant and seal it off
        if (learner$currentFit$V < V_tol) {
          learner$isLocked = TRUE
          learner$fitFunction = fitFnConstComp
          learner$predFunction = predFnConstComp
          learner$updateFit()
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
        learner_name = paste("mu_", learner$root$nextNumber, sep = '')
        combine_name = paste("combine_", learner$root$nextNumber, sep = '')
        try({ # weird error, try to fix later, but still works
          learner$root$nextNumber = learner$root$nextNumber + 1
        }, silent = T)
        
        add_fit = list(mu1 = rep(0, length(learner$Y)), mu2 = rep(0, length(learner$Y)), KL_div = 0)
        add_node = VEBBoostNodeComp$new(learner_name, fitFunction = fitFn, predFunction = predFn, currentFit = add_fit)
        learner$AddSiblingVEB(add_node, "+", combine_name)

        learner_name = paste("mu_", learner$root$nextNumber, sep = '')
        combine_name = paste("combine_", learner$root$nextNumber, sep = '')
        try({ # weird error, try to fix later, but still works
          learner$root$nextNumber = learner$root$nextNumber + 1
        }, silent = T)
        
        mult_fit = list(mu1 = rep(1, length(learner$Y)), mu2 = rep(1, length(learner$Y)), KL_div = 0)
        mult_node = VEBBoostNodeComp$new(learner_name, fitFunction = fitFn, predFunction = predFn, currentFit = mult_fit)
        learner$children[[1]]$AddSiblingVEB(mult_node, "*", combine_name)
      }

      invisible(self$root)
    },
    
    convergeFitAll = function(tol = 1e-3, V_tol = 1e-3, update_sigma2 = FALSE, update_ELBO_progress = TRUE, verbose = TRUE) {
      self$addLearnerAll(V_tol)
      self$convergeFit(tol, update_sigma2, update_ELBO_progress, verbose)
      return(invisible(self$root))
    },
    
    addLearnerEB = function(learner_name, combine_name, fitFn = NULL, predFn = NULL, tol = 1e-3, V_tol = 1e-3) {
      if (self$ebCombineOperator == ".") { # if node is "locked"
        return(invisible(self))
      } else if (self$ebCombineOperator == "+") {
        eb_current_fit = list(mu1 = rep(0, length(self$Y)), mu2 = rep(0, length(self$Y)), KL_div = 0)
      } else {
        eb_current_fit = list(mu1 = rep(1, length(self$Y)), mu2 = rep(1, length(self$Y)), KL_div = 0)
      }
      
      node_clone = VEBBoostNodeComp$new(self$name, fitFunction = self$fitFunction, predFunction = self$predFunction, currentFit = self$currentFit)
      node_clone$X = self$X
      node_clone$Y = self$Y
      node_clone$sigma2 = self$sigma2
      
      if (is.null(fitFn)) {
        fitFn = self$fitFunction
      }
      if (is.null(predFn)) {
        predFn = self$predFunction
      }
      
      eb_node = VEBBoostNodeComp$new(learner_name, fitFunction = fitFn, predFunction = predFn, currentFit = eb_current_fit)
      eb_learner = node_clone$AddSiblingVEB(eb_node, self$ebCombineOperator, "eb_combine")
      eb_learner$convergeFit(tol)
      
      if (eb_learner$children[[2]]$currentFit$V < V_tol) { # if adding basically a constant
        if (self$ebCombineOperator == "+") { # try multiplicative combine
          self$ebCombineOperator = "*"
          return(self$addLearnerEB(learner_name, combine_name, fitFn, predFn, tol, V_tol))
        } else { # add multiplicative constant, "lock" node
          eb_node$fitFunction = fitFnConstComp
          eb_node$predFunction = predFnConstComp
          eb_learner$convergeFit(tol)
          self$currentFit = eb_learner[[self$name]]$currentFit
          ebCombineOperator = self$ebCombineOperator
          self$ebCombineOperator = "."
          new_node = VEBBoostNodeComp$new(learner_name, fitFunction = eb_node$fitFunction, predFunction = eb_node$predFunction, currentFit = eb_node$currentFit, ebCombineOperator = ".")
          self$AddSiblingVEB(new_node, ebCombineOperator, combine_name)
          self$updateMoments()
        }
      } else {
        ebCombineOperator = self$ebCombineOperator
        self$ebCombineOperator = "+"
        self$currentFit = eb_learner[[self$name]]$currentFit
        new_node = VEBBoostNodeComp$new(learner_name, fitFunction = eb_node$fitFunction, predFunction = eb_node$predFunction, currentFit = eb_node$currentFit)
        self$AddSiblingVEB(new_node, ebCombineOperator, combine_name)
        self$updateMoments()
      }
      
      invisible(self)
    }, 
    
    convergeFitEB = function(tol = 1e-3, V_tol = 1e-3, update_sigma2 = FALSE, update_ELBO_progress = TRUE, verbose = TRUE) {
      base_learners = Traverse(self$root, filterFun = function(x) x$isLeaf & (x$ebCombineOperator != "."))
      for (learner in base_learners) { # change any near-constant leaf to a constant and seal it off
        if (learner$currentFit$V < V_tol) {
          learner$fitFunction = fitFnConstComp
          learner$predFunction = predFnConstComp
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
    
    # jointly add to all nodes
    # convergeFitEBJoint = function(tol = 1e-3, V_tol = 1e-3, update_sigma2 = FALSE, update_ELBO_progress = TRUE, verbose = TRUE) {
    #   base_learners = Traverse(self$root, filterFun = function(x) x$isLeaf & (x$ebCombineOperator != "."))
    #   for (learner in base_learners) { # change any near-constant leaf to a constant and seal it off
    #     if (learner$currentFit$V < V_tol) {
    #       learner$fitFunction = fitFnConstComp
    #       learner$ebCombineOperator = "."
    #       learner$updateFit()
    #     }
    #   }
    #   
    #   base_learners = Traverse(self$root, filterFun = function(x) x$isLeaf & (x$ebCombineOperator != "."))
    #   for (learner in base_learners) { # add to all nodes
    #     learner_name = paste("mu_", learner$root$nextNumber, sep = '')
    #     combine_name = paste("combine_", learner$root$nextNumber, sep = '')
    #     learner$addLearnerEB(learner_name, combine_name)
    #     try({ # weird error, try to fix later, but still works
    #       learner$root$nextNumber = learner$root$nextNumber + 1
    #     }, silent = T)
    #   }
    #   
    #   # update fit
    #   ELBOs = numeric(1000)
    #   ELBOs[1] = -Inf
    #   ELBOs[2] = self$root$ELBO
    #   i = 2
    #   while (abs(ELBOs[i] - ELBOs[i-1]) > tol) {
    #     self$root$Do(function(x) x$updateFit(), filterFun = function(x) x$isLeaf)
    #     if (update_sigma2) {
    #       self$root$update_sigma2()
    #     }
    #     i = i + 1
    #     if (i > length(ELBOs)) { # double size of ELBOs for efficiency rather than making it bigger each iteration
    #       ELBOs = c(ELBOs, rep(0, i))
    #     }
    #     ELBOs[i] = self$root$ELBO
    #     if (verbose & ((i %% 100) == 0)) {
    #       cat(paste("ELBO: ", ELBOs[i], sep = ""))
    #       cat("\n")
    #     }
    #   }
    #   ELBOs = ELBOs[2:i]
    #   
    #   if (update_ELBO_progress) {
    #     self$ELBO_progress = ELBOs
    #   }
    #   
    #   invisible(self$root)
    # }, 
    
    # predict.veb = function(X_new) { # function to get prediction on new data
    #   self$root$Do(function(node) {
    #     node$pred_mu1 = node$predFunction(X_new, node$currentFit, 1)
    #     node$pred_mu2 = node$predFunction(X_new, node$currentFit, 2)
    #   }, filterFun = function(x) x$isLeaf)
    #   invisible(self$root)
    # }
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
    .mu1 = NA, # current first moment
    .mu2 = NA, # current second moment
    .ELBO_progress = list(-Inf), # list of ELBOs for each iteration of growing and fitting the tree
    .pred_mu1 = NULL, # prediction based on predFunction and given new data (first moment)
    .pred_mu2 = NULL, # prediction based on predFunction and given new data (second moment)
    .isLocked = FALSE # locked <=> V < V_tol, or both learners directly connected (sibling and parent's sibling) are constant
  ), 
  
  active = list(
    
    mu1 = function(value) {
      if (!missing(value)) {
        stop("`$mu1` cannot be modified directly", call. = FALSE)
      }
      if (self$isLeaf) {
        return(self$currentFit$mu1)
      } else {
        return(private$.mu1)
      }
    },
    
    mu2 = function(value) {
      if (!missing(value)) {
        stop("`$mu2` cannot be modified directly", call. = FALSE)
      }
      if (self$isLeaf) {
        return(self$currentFit$mu2)
      } else {
        return(private$.mu2)
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
          if (self$root$family == "gaussian") {
            return(private$.Y)
          } else if (self$root$family == "binomial") {
            g_xi = g(self$root$xi) # vector of g(xi_i), pre-compute once
            d = ((g_xi - .5) / self$root$xi) # vector of (g(xi_i) - .5) / xi_i, pre-compute once
            d[self$root$xi == 0] = .25 # case of 0/0 (e.g. x_i is all 0), use L'Hopital
            return((private$.Y - .5) / d)
          } else {
            stop("family must be either 'gaussian' or 'binomial")
          }
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
            g_xi = g(self$root$xi) # vector of g(xi_i), pre-compute once
            d = ((g_xi - .5) / self$root$xi) # vector of (g(xi_i) - .5) / xi_i, pre-compute once
            d[self$root$xi == 0] = .25 # case of 0/0 (e.g. x_i is all 0), use L'Hopital
            return(1 / d)
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
        g_xi = g(self$root$xi) # vector of g(xi_i), pre-compute once
        d = ((g_xi - .5) / self$root$xi) # vector of (g(xi_i) - .5) / xi_i, pre-compute once
        d[self$root$xi == 0] = .25 # case of 0/0 (e.g. x_i is all 0), use L'Hopital
        return(
          sum(log(g_xi)) + sum((self$root$xi / 2)*(d*self$root$xi - 1)) + sum((self$root$raw_Y - .5) * self$root$mu1) - 
            .5*sum(self$root$mu2 * d) - self$KL_div
        )
      } else {
        stop("family must be either 'gaussian' or 'binomial")
      }
    }, 

    xi = function(value) { # optimal variational parameters, set to +sqrt(mu2)
      if (!missing(value)) {
        stop("`$xi` cannot be modified directly", call. = FALSE)
      }
      return(sqrt(self$root$mu2))
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
fitFnConstComp = function(X, Y, sigma2, init) {
  if (length(sigma2) == 1) {
    sigma2 = rep(sigma2, length(Y))
  }
  intercept = weighted.mean(Y, 1/sigma2)
  KL_div = 0

  mu1 = rep(intercept, length(Y))
  mu2 = rep(intercept^2, length(Y))
  return(list(mu1 = mu1, mu2 = mu2, intercept = intercept, KL_div = KL_div))
}

predFnConstComp = function(X_new, currentFit, moment = c(1, 2)) {
  if (moment == 1) {
    return(rep(currentFit$intercept, length(currentFit$mu1)))
  } else if (moment == 2) {
    return(rep(currentFit$intercept^2, length(currentFit$mu2)))
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
  
