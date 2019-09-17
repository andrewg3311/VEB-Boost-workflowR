### Functions for fitting a given VEB-Tree ###
#' Each node in the tree stores the following values:
#' 1. First and second moments of the current fit of the model defined by the subtree w/ self as the root
#' 2. Sum of the KL divergence of all base learners that are descendants
#' 
#' Internal nodes also store the operator at that node, either "+" or "*"
#' 
#' Terminal nodes also must store the following:
#' 1. A function that solves the weighted VEB regression problem (can initialize from previous fit, "warm start")
#' 2. A current fit, which has fields $mu (for 1st moment), $mu2 (for second moment), and $KL_div (for KL divergence)
#' 
#' Root also stores:
#' 1. Response Y
#' 2. Current value for sigma^2


library(data.tree)

# function to get a numeric vector for the path to the given node
getPathVector = function(node, path_vector = numeric()) {
  if (node$isRoot) {
    return(path_vector)
  }
  return(c(getPathVector(node$parent, path_vector), node$position, path_vector))
}

# function to get the pseudo data used to update the node  at path_vectir (response and variances)
# subtree is sub-tree
getPseudoData = function(subtree, path_vector, Y_tilde = node$root$Y, sigma2_tilde = node$root$sigma2) {
  
  if (length(path_vector) == 0) { # if for some reason only one learner in tree, needed to avoid potential subscript out-of-bounds errors
    return(list(Y_tilde = Y_tilde, sigma2_tilde = sigma2_tilde))
  }
  
  other_half = subtree$children[[3 - path_vector[1]]] # requires binary tree, half excluding base learner
  same_half = subtree$children[[path_vector[1]]] # half including desired base learner
  
  if (subtree$operator == "+") { # add adding left and right halves
    Y_tilde = Y_tilde - other_half$mu
  } else { # if Schur multiplying left and right halves
    Y_tilde = Y_tilde * (other_half$mu / other_half$mu2)
    sigma2_tilde = sigma2_tilde / other_half$mu2
  }
  
  if (length(path_vector) == 1) { # if at parent of desired base learner node
   return(list(Y_tilde = Y_tilde, sigma2_tilde = sigma2_tilde))
  }
  
  return(getPseudoData(subtree$children[[path_vector[1]]], path_vector[-1], Y_tilde, sigma2_tilde))
  
}

# function to update the base learner at the terminal node
# node has a function stored in $fit_function
# takes in data (and initialize point)
# stored in $current_fit
updateBaseLearner = function(terminalNode) {
  
  path_vector = getPathVector(terminalNode) # get numeric path vector
  pseudo_data = getPseudoData(terminalNode$root, path_vector) # get pseudo data
  
  terminalNode$current_fit = fit_function(Y = pseudo_data$Y_tilde, sigma2 = pseudo_data$sigma2_tilde,
                                          init = ifelse(is.null(terminalNode$current_fit), NULL, terminalNode$current_fit))
  
  # update stored values for moments and KL divergence
  terminalNode$mu = terminalNode$current_fit$mu
  terminalNode$mu2 = terminalNode$current_fit$mu2
  terminalNode$KL_div = terminalNode$current_fit$KL_div
  
  if (terminalNode$isRoot) { # if only 1 learner in tree
    return(terminalNode)
  }
  
  # pass updated values up the tree
  parent_node = terminalNode$parent
  while (length(path_vector) > 0) {
    
    other_half = parent_node$children[[3 - path_vector[length(path_vector)]]] # requires binary tree, other half of tree
    same_half = parent_node$children[[path_vector[1]]] # half including desired base learner
    
    if (parent_node$operator == "+") { # if adding sub-learners
      parent_node$mu = same_half$mu + other_half$mu
      parent_node$mu2 = same_half$mu2 + other_half$mu2
    } else { # if multiplying sub-learners
      parent_node$mu = same_half$mu * other_half$mu
      parent_node$mu2 = same_half$mu2 * other_half$mu2
    }
    
    parent_node$KL_div = same_half$KL_div + other_half$KL_div # add KL divergences
    
    path_vector = path_vector[-length(path_vector)]
    parent_node = parent_node$parent
  }
  
  return(terminalNode)
  
}




# function to update all base learners
updateAll = function(tree) {
  tree$Do(updateBaseLearner, filterFun = function(x) x$isLeaf)
}
