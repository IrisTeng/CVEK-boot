gpr <- 
  function(formula, label.names, 
           data = NULL, method = NULL, l = 1, 
           lambda = exp(seq(-5, 5, 1))){
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
    #   null.hypo: logical. To indicate whether the Kernel function is calculated
    #              under the null hypothesis.
    #
    # Returns:
    #   sigma.n: A numeric number chosen as the best size of noise via LooCV.
    #   intercept: A numberic number indicating the bias weight.
    #   alpha: A numeric string indicating the projection vector.
    
    # prepare data
    y <- data[, as.character(attr(terms(formula), "variables"))[2]]
    
    # extract the information from the given formula to construct a true formula
    re <- GenericFormula(formula, label.names)
    generic.formula0 <- re[[1]]
    len <- re[[2]]
    
    # extract design matrix
    X <- model.matrix(generic.formula0, data)[, -1]
    n <- nrow(X)
    
    # estimation
    Kern <- KernelGenerate(method, l)
    X1 <- X[, c(1:length(label.names[[1]]))]
    X2 <- X[, c((length(label.names[[1]]) + 1):len)]
    K <- Kern(X1, X1) + Kern(X2, X2)
    
    if (length(lambda) == 1)
      lambda0 = 0
    else 
      lambda0 <- LooCV(y, K, lambda)[[1]]

    K1 <- cbind(1, K)
    K2 <- cbind(0, rbind(0, K))
    
    theta <- ginv(lambda0 * K2 + t(K1) %*% K1) %*% t(K1) %*% y
    beta0 <- theta[1]
    alpha <- theta[-1]

    return(list(sigma2.n = lambda0, intercept = beta0, alpha = alpha, K = K))
  }
