### Node Object ###

VEBBoostNode <- R6::R6Class(
  "VEBBoostNode", 
  public = list(

  ),
  active = list(
    mu1 = 0, # first moment of learner defined by sub-tree with this node as the root
    
    mu2 = 0, # second moment of learner defined by sub-tree with this node as the root
    #position = NA # gives the position in the tree
    # NULL <=> isRoot
    # represented as a bit-string (logical vector)
    # 0/FALSE means left, 1/TRUE means right
    
    KL_div = 0, # KL divergence from q to g
    
    Y = function() { # response for sub-tree
      if (self$isRoot) {
        return(self$Y)
      }
      if (self$parent$operator == "+") {
        return(self$parent$Y - self$siblings[[1]]$mu1)
      }
      if (self$parent$operator == "*") {
        return(self$parent$Y * (self$siblings[[1]]$mu1 / self$siblings[[1]]$mu2))
      }
    },
    
    sigma2 = function() { # variance for sub-tree
      if (self$isRoot) {
        return(self$sigma2)
      }
      if (self$parent$operator == "+") {
        return(self$parent$sigma2)
      }
      if (self$parent$operator == "*") {
        return(self$parent$sigma2 / self$siblings[[1]]$mu2)
      }
    }
  ), 
  inherit = data.tree::Node, 
  lock_objects = FALSE
)
