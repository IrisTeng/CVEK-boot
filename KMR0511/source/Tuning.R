Tuning <- 
  function(y, K.mat, lambda = exp(seq(-5, 5, 1)), mode = 'loocv'){
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
    if (mode == "loocv"){
      CV <- sapply(lambda, function(k){
        A <- K.mat %*% ginv(K.mat + k * diag(n))
        sum(((diag(n) - A) %*% y / diag(diag(n) - A)) ^ 2)
      })
    }
    else if (mode == "AICc"){
      CV <- sapply(lambda, function(k){
        A <- K.mat %*% ginv(K.mat + k * diag(n))
        log(t(y) %*% (diag(n) - A) %*% (diag(n) - A) %*% y) +
          2 * (tr(A) + 2) / (n - tr(A) - 3)
      })
    }
    else if (mode == "GCVc"){
      CV <- sapply(lambda, function(k){
        A <- K.mat %*% ginv(K.mat + k * diag(n))
        log(t(y) %*% (diag(n) - A) %*% (diag(n) - A) %*% y) -
          2 * log(max(0, 1 - tr(A) / n - 2 / n))
      })
    }
    else if (mode == "gmpml"){
      CV <- sapply(lambda, function(k){
        A <- K.mat %*% ginv(K.mat + k * diag(n))
        log(t(y) %*% (diag(n) - A) %*% y) - 
          1 / (n - 1) * log(det((diag(n) - A)))
      })
    }
    else
      stop("mode must be loocv, AICc, GCVc or gmpml!")
    
    lambda0 <- lambda[which(CV == min(CV))]
    # A <- K.mat %*% ginv(K.mat + lambda0 * diag(n))
    # error <- (diag(n) - A) %*% y / diag(diag(n) - A)
    return(lambda0)
  }
