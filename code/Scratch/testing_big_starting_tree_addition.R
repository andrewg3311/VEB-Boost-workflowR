learner_add = VEBBoostNodeComp$new("mu_0", fitFunction = fitFn2_stumps_int, predFunction = predFn2_stumps, currentFit = list(mu1 = rep(0, length(Y)), mu2 = rep(0, length(Y)), KL_div = 0))
learner_add$X = X
learner_add$Y = Y
learner_add$sigma2 = 1

for (i in 1:8) {
  learners = Traverse(learner_add, filterFun = function(x) x$isLeaf)
  for (learner in learners) {
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
  }
}
 
ELBOs = numeric(1000)
ELBOs[1] = learner_add$ELBO

base_learners = Traverse(learner_add, filterFun = function(x) x$isLeaf)

i = 2
system.time({
  for (learner in base_learners) {
    learner$updateFit()
  
    ELBOs[i] = learner_add$ELBO
    i = i + 1
    if (i == 62) {
      break
    }
  }
})
learner_add$update_sigma2()
ELBOs[i] = learner_add$ELBO

i = i + 1
system.time({
  for (learner in base_learners) {
    learner$updateFit()
    
    ELBOs[i] = learner_add$ELBO
    i = i + 1
  }
})
learner_add$update_sigma2()
ELBOs[i] = learner_add$ELBO

plot(ELBOs[1:i])



learner_add2 = VEBBoostNodeComp$new("mu_0", fitFunction = fitFn2_stumps, predFunction = predFn2_stumps, currentFit = list(mu1 = rep(0, length(Y)), mu2 = rep(0, length(Y)), KL_div = 0))
learner_add2$X = X
learner_add2$Y = Y
learner_add2$sigma2 = 1

for (i in 1:8) {
  learners = Traverse(learner_add2, filterFun = function(x) x$isLeaf)
  for (learner in learners) {
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
  }
}

ELBOs2 = numeric(1000)
ELBOs2[1] = learner_add2$ELBO

base_learners2 = Traverse(learner_add2, filterFun = function(x) x$isLeaf)

i = 2
system.time({
  for (learner in base_learners2) {
    learner$updateFit()
    
    ELBOs2[i] = learner_add2$ELBO
    i = i + 1
  }
})
learner_add2$update_sigma2()
ELBOs2[i] = learner_add2$ELBO

i = i + 1
system.time({
  for (learner in base_learners2) {
    learner$updateFit()
    
    ELBOs2[i] = learner_add2$ELBO
    i = i + 1
  }
})
learner_add2$update_sigma2()
ELBOs2[i] = learner_add2$ELBO

plot(ELBOs2[1:i])


learner_add3 = VEBBoostNodeComp$new("mu_0", fitFunction = fitFn2, predFunction = predFn2, currentFit = list(mu1 = rep(0, length(Y)), mu2 = rep(0, length(Y)), KL_div = 0))
learner_add3$X = X
learner_add3$Y = Y
learner_add3$sigma2 = 1

for (i in 1:8) {
  learners = Traverse(learner_add3, filterFun = function(x) x$isLeaf)
  for (learner in learners) {
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
  }
}

ELBOs3 = numeric(1000)
ELBOs3[1] = learner_add3$ELBO

base_learners3 = Traverse(learner_add3, filterFun = function(x) x$isLeaf)

i = 2
system.time({
  for (learner in base_learners3) {
    learner$updateFit()
    
    ELBOs3[i] = learner_add3$ELBO
    i = i + 1
  }
})
learner_add3$update_sigma2()
ELBOs3[i] = learner_add3$ELBO

i = i + 1
system.time({
  for (learner in base_learners3) {
    learner$updateFit()
    
    ELBOs3[i] = learner_add3$ELBO
    i = i + 1
  }
})
learner_add3$update_sigma2()
ELBOs3[i] = learner_add3$ELBO

plot(ELBOs3[1:i])



learner_add4 = VEBBoostNodeComp$new("mu_0", fitFunction = fitFn2_int, predFunction = predFn2, currentFit = list(mu1 = rep(0, length(Y)), mu2 = rep(0, length(Y)), KL_div = 0))
learner_add4$X = X
learner_add4$Y = Y
learner_add4$sigma2 = 1

for (i in 1:8) {
  learners = Traverse(learner_add4, filterFun = function(x) x$isLeaf)
  for (learner in learners) {
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
  }
}

ELBOs4 = numeric(1000)
ELBOs4[1] = learner_add4$ELBO

base_learners4 = Traverse(learner_add4, filterFun = function(x) x$isLeaf)

i = 2
system.time({
  for (learner in base_learners4) {
    learner$updateFit()
    
    ELBOs4[i] = learner_add4$ELBO
    i = i + 1
  }
})
learner_add4$update_sigma2()
ELBOs4[i] = learner_add4$ELBO

i = i + 1
system.time({
  for (learner in base_learners4) {
    learner$updateFit()
    
    ELBOs4[i] = learner_add4$ELBO
    i = i + 1
  }
})
learner_add4$update_sigma2()
ELBOs4[i] = learner_add4$ELBO

plot(ELBOs4[1:i])



learner_add5 = VEBBoostNodeComp$new("mu_0", fitFunction = fitFnFull, predFunction = predFnFull, currentFit = list(mu1 = rep(0, length(Y)), mu2 = rep(0, length(Y)), KL_div = 0))
learner_add5$X = X
learner_add5$Y = Y
learner_add5$sigma2 = 1

for (i in 1:8) {
  learners = Traverse(learner_add5, filterFun = function(x) x$isLeaf)
  for (learner in learners) {
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
  }
}

ELBOs5 = numeric(1000)
ELBOs5[1] = learner_add5$ELBO

base_learners5 = Traverse(learner_add5, filterFun = function(x) x$isLeaf)

i = 2
system.time({
  for (learner in base_learners5) {
    learner$updateFit()
    
    ELBOs5[i] = learner_add5$ELBO
    i = i + 1
  }
})
learner_add5$update_sigma2()
ELBOs5[i] = learner_add5$ELBO

i = i + 1
system.time({
  for (learner in base_learners5) {
    learner$updateFit()
    
    ELBOs5[i] = learner_add5$ELBO
    i = i + 1
  }
})
learner_add5$update_sigma2()
ELBOs5[i] = learner_add5$ELBO

plot(ELBOs5[1:i])


