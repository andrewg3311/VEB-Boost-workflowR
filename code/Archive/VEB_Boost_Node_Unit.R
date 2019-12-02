### Unit Node Object ###

VEBBoostNodeUnit <- R6::R6Class(
  "VEBBoostNodeUnit", 
  public = list(
    fitFunction = function(Y, sigma2, init) {
      return(mu1 = rep(1, length(Y)), mu2 = rep(1, length(Y)), KL_div = 0)
    }
    
  ),
  inherit = VEBBoostNode, 
  lock_objects = FALSE
)
