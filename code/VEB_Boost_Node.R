### Node Object ###

VEBBoostNode <- R6::R6Class(
  "VEBBoostNode", 
  public = list(
    operator = NULL, # either "+" or "*" for internal nodes, NULL for terminal nodes
    
    fitFunction = NULL, # function that takes in response Y, variances sigma2, and returns the fit
    # the fit must have fields mu1 (first moment), mu2 (second moment), and KL_div (KL divergence from q to g)
    
    currentFit = NULL, # current fit for fitting function
    
    # initialize = function(..., operator = NULL, fitFunction = NULL, currentFit = NULL) {
    #   self = data.tree::Node$new(...)
    #   self$operator = operator
    #   self$fitFunction = fitFunction
    #   self$currentFit = currentFit
    #   
    #   # checks
    #   if (self$isLeaf) { # if terminal node
    #     if (!is.null(self$operator)) {
    #       stop("Terminal nodes cannot have an operator")
    #     }
    #     if (is.null(fitFunction)) {
    #       stop("Terminal nodes must have a fit function")
    #     }
    #   } else { # if internal node
    #     if (is.null(self$operator)) {
    #       stop("Internal nodes must have an operator")
    #     } else if (!(operator %in% c("+", "*"))) {
    #       stop("Invalid operator, must be one of '+' or '*'")
    #     }
    #     if (!is.null(self$fitFunction)) {
    #       stop("Internal nodes cannot have a fit function")
    #     }
    #   }
    # 
    # },
    
    AddChildVEB = function(name, check = c("check", "no-warn", "no-check"), ...) { # add VEB node as child
      child = VEBBoostNode$new(as.character(name), check, ...)
      invisible(self$AddChildNode(child))
    },
    
    updateFit = function() { # function to update currentFit
      self$currentFit = self$fitFunction(Y = self$Y, sigma2 = self$sigma2, init = self$currentFit)
      invisible(self)
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
    }
  ), 
  inherit = data.tree::Node, 
  lock_objects = FALSE
)
