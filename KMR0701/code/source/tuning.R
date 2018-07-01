tuning <- 
  function(Y, K_mat, mode, lambda){
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
    n <- nrow(K_mat)
    
    # estimation
    if (mode == "loocv"){
      CV <- sapply(lambda, function(k){
        A <- K_mat %*% ginv(K_mat + k * diag(n))
        sum(((diag(n) - A) %*% Y / diag(diag(n) - A)) ^ 2)
      })
    }
    else if (mode == "AICc"){
      CV <- sapply(lambda, function(k){
        A <- K_mat %*% ginv(K_mat + k * diag(n))
        log(t(Y) %*% (diag(n) - A) %*% (diag(n) - A) %*% Y) +
          2 * (tr(A) + 2) / (n - tr(A) - 3)
      })
    }
    else if (mode == "GCVc"){
      CV <- sapply(lambda, function(k){
        A <- K_mat %*% ginv(K_mat + k * diag(n))
        log(t(Y) %*% (diag(n) - A) %*% (diag(n) - A) %*% Y) -
          2 * log(max(0, 1 - tr(A) / n - 2 / n))
      })
    }
    else if (mode == "gmpml"){
      CV <- sapply(lambda, function(k){
        A <- K_mat %*% ginv(K_mat + k * diag(n))
        log(t(Y) %*% (diag(n) - A) %*% Y) - 
          1 / (n - 1) * log(det((diag(n) - A)))
      })
    }
    else
      stop("mode must be loocv, AICc, GCVc or gmpml!")
    
    lambda0 <- lambda[which(CV == min(CV))]
    return(lambda0)
  }
