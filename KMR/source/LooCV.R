LooCV <- 
  function(formula, label.names, 
           data = NULL, method = NULL, l = 1, lambda = exp(seq(-5, 5, 1))){
    # An implementation of Gaussian processes for regression.
    #
    # Args:
    #   formula:  A symbolic description of the model to be fitted.
    #   label.names: A character string indicating all the interior variables
    #                included in each predictor.
    #   data: A dataframe to be fitted.
    #   method: A character string indicating which Kernel is to be computed.
    #   l: A numeric number indicating the hyperparameter of a specific kernel.
    #   lambda: A numeric string specifying the range of noise to be chosen.
    #           The lower limit of lambda must be above 0.
    #
    # Returns:
    #   lambda0: A numeric number chosen as the best size of noise via LooCV.
    
    # prepare data
    y <- data[, as.character(attr(terms(formula), "variables"))[2]]
    
    # extract the information from the given formula to construct a true formula
    generic.formula0 <- 
      GenericFormula(formula, label.names)$generic.formula
    
    # extract design matrix
    X <- model.matrix(generic.formula0, data=data)[, -1]
    n <- nrow(X)
    
    # estimation
    Kern <- KernelGenerate(method, l)
    K <- Kern(X, X)
    CV <- sapply(lambda, function(k){
      A <- K %*% ginv(K + k * diag(n))
      sum(((diag(n) - A) %*% y / diag(diag(n) - A)) ^ 2)
    })
    lambda0 <- lambda[which(CV == min(CV))]
    return(lambda0)
  }