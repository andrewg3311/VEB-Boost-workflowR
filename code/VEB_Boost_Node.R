### Node Object ###

VEBBoostNode <- R6::R6Class(
  "VEBBoostNode", 
  public = list(
    operator = NULL, # either "+" or "*" for internal nodes, NULL for terminal nodes
    
    fitFunction = NULL, # function that takes in response Y, variances sigma2, and returns the fit
    # the fit must have fields mu1 (first moment), mu2 (second moment), and KL_div (KL divergence from q to g)
    
    currentFit = NULL, # current fit for fitting function
    
    #ELBO_progress = list(), # track ELBOs for each tree fit along the way
    
    AddChildVEB = function(name, check = c("check", "no-warn", "no-check"), ...) { # add VEB node as child
      child = VEBBoostNode$new(as.character(name), check, ...)
      invisible(self$AddChildNode(child))
    },
    
    updateFit = function() { # function to update currentFit
      self$currentFit = self$fitFunction(Y = self$Y, sigma2 = self$sigma2, init = self$currentFit)
      invisible(self)
    }, 
    
    convergeFit = function(tol = 1e-6) {
      ELBOs = numeric(100)
      ELBOs[1] = -Inf
      ELBOs[2] = self$root$ELBO
      i = 2
      while (ELBOs[i] - ELBOs[i-1] > tol) {
        self$root$Do(function(x) x$updateFit(), filterFun = function(x) x$isLeaf)
        i = i+1
        ELBOs[i] = self$root$ELBO
      }
      ELBOs = ELBOs[2:i]
      
      #self$root$ELBO_progress[[length(self$root$ELBO_progress) + 1]] = ELBOs
      invisible(self$root)
    }, 
    
    AddSiblingVEB = function(learner, operator = c("+", "*"), name) { # function to add subtree as sibling to given node, combining with operator
      # learner is tree to add to self
      # operator is how to combine them
      # name is what to call the combining node
      
      # if (self$isRoot) { # if adding at root
      #   newRoot = VEBBoostNode$new(name, operator = "operator")
      #   newRoot$AddChildNode(self)
      #   newRoot$AddChildNode(learner)
      # } else {
      #   # add combining node
      #   self$parent$AddChildVEB(name, operator = "+")
      #   
      #   # fix child/parent pointers for self
      #   self$parent$children$name$private$p_children[[self$name]] = self
      #   self$parent$children$name[[self$name]] = self
      #   self$parent = self$parent$children$name
      #   
      #   # add new learner
      #   self$parent$AddChildNode(learner)
      # }
      # 
      # invisible(self$parent)
      
      learner_copy = self$clone()
      
      self$fitFunction = NULL
      self$currentFit = NULL
      self$operator = operator
      self$name = name
      
      self$AddChildNode(learner_copy)
      self$AddChildNode(learner)
      
      # combineNode = VEBBoostNode$new(name, operator = operator) # new combining internal node
      # if (!self$isRoot) { # add new node to existing tree if needed
      #   self$parent$AddChildNode(combineNode)
      # }
      # # make new node the parent of self and learner to combine with
      # 
      # if (self$isRoot) { # remove private Y and sigma2 if this was the root (shouldn't matter, just a precaution)
      #   combineNode$Y = self$Y
      #   combineNode$sigma2 = self$sigma2
      #   #combineNode$ELBO_progress = self$ELBO_progress
      #   self$Y = NA
      #   self$sigma2 = NA
      #   #self$ELBO_progress = list()
      # } else {
      #   self$parent$p_children[[self$name]] = NULL
      #   self$parent[[self$name]] = NULL
      #   self$parent = NULL
      # }
      # if (learner$isRoot) {
      #   learner$Y = NA
      #   learner$sigma2 = NA
      # } else {
      #   learner$parent$p_children[[learner$name]] = NULL
      #   learner$parent[[learner$name]] = NULL
      #   learner$parent = NULL
      # }
      # combineNode$AddChildNode(self)
      # combineNode$AddChildNode(learner)
      
      invisible(self$root) # return root, so we can assign entire tree to new value (in case we're splitting at root)
      
    }, 
    
    addLearner = function(learner_name, combine_name, fitFn = NULL, currentFit = NULL, tol = 1e-3) { # add a new learner to the tree
      base_learners = Traverse(self$root, filterFun = function(x) x$isLeaf)
      var_inf = sapply(base_learners, function(x) x$variance_influence) # vector of var inf values
      node_to_split = base_learners[[which.max(var_inf)]] # note: if all additive, this will lead to an unbalanced tree (so calculating moments, etc not as efficient)
      
      node_clone = VEBBoostNode$new(node_to_split$name, fitFunction = node_to_split$fitFunction, currentFit = node_to_split$currentFit)
      node_clone$Y = node_to_split$Y
      node_clone$sigma2 = node_to_split$sigma2
      node_clone2 = node_clone$clone()
      
      if (is.null(fitFn)) {
        fitFn = node_clone$fitFunction
      }
      if (is.null(currentFit)) {
        currentFitAdd = list(mu1 = rep(0, length(node_clone$Y)), mu2 = rep(1, length(node_clone$Y)), KL_div = 0)
        currentFitMult = list(mu1 = rep(1, length(node_clone$Y)), mu2 = rep(1, length(node_clone$Y)), KL_div = 0)
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
    }
    
    
  ),
  private = list(
    .Y = NA, # response, only not NA at root
    .sigma2 = NA # variance, only not NA at root
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
      } else if (self$operator == "*") {
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
      } else if (self$operator == "*") {
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
    
    ELBO = function(value) { # ELBO for entire tree
      if (!missing(value)) {
        stop("`$ELBO` cannot be modified directly", call. = FALSE)
      }
      # make s2 vector of variances
      s2 = ifelse(length(self$root$sigma2) == 1, rep(self$root$sigma2, length(self$root$Y)), self$root$sigma2)
      return(
        (-.5 * sum(log(2*pi*s2))) - 
          (.5 * (sum(self$root$Y^2 / s2) - 2*sum(self$root$Y * self$root$mu1 / s2) + sum(self$root$mu2 / s2))) - 
          self$root$KL_div
      )
    }, 
    
    variance = function(value) { # variance of mean values at node
      if (!missing(value)) {
        stop("`$variance` cannot be modified directly", call. = FALSE)
      }
      
      return(self$mu2 - self$mu1^2)
    }, 
    
    variance_partial_deriv = function(value) { # partial derivative of variance of mean w.r.t. variance at node
      if (!missing(value)) {
        stop("`$variance_partial_deriv` cannot be modified directly", call. = FALSE)
      }
      
      if (self$isRoot) {
        return(rep(1, length(self$Y)))
      }
      if (self$parent$operator == "+") {
        return(self$parent$variance_partial_deriv)
      } else {
        return(self$siblings[[1]]$variance * self$parent$variance_partial_deriv)
      }
    }, 
    
    variance_influence = function(value) { # weighted variance partial deriv, weighted by variance at node
      if (!missing(value)) {
        stop("`$variance_influence` cannot be modified directly", call. = FALSE)
      }
      
      return(weighted.mean(self$variance, self$variance_partial_deriv))
    }
  ), 
  inherit = data.tree::Node, 
  lock_objects = FALSE
)
