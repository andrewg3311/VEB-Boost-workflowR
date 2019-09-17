### Terminal Node Object ###

VEBTerminalNode <- R6::R6Class(
  "VEBTerminalNode", 
  public = list(
    
  ),
  private = list(
    fitFunction = NULL, # function that takes in predictors X, response Y, variances sigma2, and returns the fit
    # the fit must have fields mu1 (first moment), mu2 (second moment), and KL_div (KL divergence from q to g)
    currentFit = NULL # current fit for fitting function
  ), 
  active = list(
    mu1 = function() {
      self$currentFit$mu1
    },
    
    mu2 = function() {
      self$currentFit$mu2
    },
    
    KL_div = function() {
      self$currentFit$KL_div
    }
  ), 
  inherit = VEBBoostNode, 
  lock_objects = FALSE
)
