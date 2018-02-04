LooCV <- 
  function(y, K.mat, lambda = exp(seq(-5, 5, 1))){
    # An implementation of Gaussian processes for regression.
    #
    # Args:
    #   y: A vector of response from original data.
    #   K.mat: Kernel matrix calculated from the original data.
    #   lambda: A numeric string specifying the range of noise to be chosen.
    #           The lower limit of lambda must be above 0.
    #
    # Returns:
    #   lambda0: A numeric number chosen as the best size of noise via LooCV.
    
    # prepare data
    n <- nrow(K.mat)
    # estimation
    CV <- sapply(lambda, function(k){
      A <- K.mat %*% ginv(K.mat + k * diag(n))
      sum(((diag(n) - A) %*% y / diag(diag(n) - A)) ^ 2)
    })
    lambda0 <- lambda[which(CV == min(CV))]
    
    return(lambda0)
  }
