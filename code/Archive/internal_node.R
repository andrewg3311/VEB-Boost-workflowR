### Internal Node Object ###

VEBInternalNode <- R6::R6Class(
  "VEBInternalNode", 
  public = list(
    operator = NULL, # either "+" or "*"
    # updateMoments = function() { # update moments from children
    #   if (self$operator == "+") { # if adding sub-models
    #     self$mu1 = apply(sapply(self$children, function(x) x$mu1), MARGIN = 1, sum)
    #     self$mu2 = apply(sapply(self$children, function(x) x$mu2), MARGIN = 1, sum)
    #   } else if (self$operator == "*") { # if (Schur) multiplying sub-modls
    #     self$mu1 = apply(sapply(self$children, function(x) x$mu1), MARGIN = 1, prod)
    #     self$mu2 = apply(sapply(self$children, function(x) x$mu2), MARGIN = 1, prod)
    #   } else {
    #     stop("Internal Node has invalid operator, must be either '+' or '*'")
    #   }
    #   return(invisible(self)) # per the Hadley R6 convention
    # },
    # updateKL = function() { # update KL divergence from children, always add together
    #   self$KL_div = sum(sapply(self$children, function(x) x$KL_div))
    #   return(invisible(self))
    # },
    AddChildVEB = function(name, check = c("check", "no-warn", "no-check"), internal = TRUE, ...) { # add VEB node as child
      if (internal) { # add internal node as child
        child = VEBInternalNode$new(as.character(name), check, ...)
      } else { # add terminal node as child
        child = VEBTerminalNode$new(as.character(name), check, ...)
      }
      invisible(self$AddChildNode(child))
    }
  ), 
  active = list(
    mu1 = function() {
      if (self$operator == "+") {
        return(apply(sapply(self$children, function(x) x$mu1), MARGIN = 1, sum))
      } else if (self$operator == "*") {
        return(apply(sapply(self$children, function(x) x$mu1), MARGIN = 1, prod))
      }
    },
    
    mu2 = function() {
      if (self$operator == "+") {
        return(apply(sapply(self$children, function(x) x$mu2), MARGIN = 1, sum))
      } else if (self$operator == "*") {
        return(apply(sapply(self$children, function(x) x$mu2), MARGIN = 1, prod))
      }
    },
    
    KL_div = function() {
      return(sum(sapply(self$children, function(x) x$KL_div)))
    }
  ), 
  inherit = VEBBoostNode, 
  lock_objects = FALSE
)
