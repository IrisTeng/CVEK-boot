ensemble <- 
  function(formula, label_names, Kernlist, data, 
           mode, strategy, beta, lambda){
    
    y <- data[, as.character(attr(terms(formula), "variables"))[2]]
    
    # extract the information from the given formula to construct a true formula
    re <- genericFormula(formula, label_names)
    generic_formula0 <- re[[1]]
    len <- re[[2]]
    
    # extract design matrix
    X <- model.matrix(generic_formula0, data)[, -1]
    n <- nrow(X)
    
    # estimation
    # D = 10
    X1 <- X[, c(1:length(label_names[[1]]))]
    X2 <- X[, c((length(label_names[[1]]) + 1):len)]
    
    D <- length(Kernlist)
    out <- baseEstimate(n, D, y, X1, X2, Kernlist, mode, lambda)
    A_hat <- out$A_hat
    error_mat <- out$error_mat
    
    if (strategy == "erm"){
      
      A <- error_mat
      B <- rep(0, n)
      E <- rep(1, D)
      F <- 1
      G <- diag(D)
      H <- rep(0, D)
      u_hat <- lsei(A, B, E = E, F = F, G = G, H = H)$X
      A_est <- u_hat[1] * A_hat[[1]]
      if(D != 1)
        for(d in 2:D)
          A_est <- A_est + u_hat[d] * A_hat[[d]]
    }
    else if (strategy == "average"){
      
	  u_hat <- rep(1 / D, D)
      A_est <- (1 / D) * A_hat[[1]]
      if(D != 1)
        for(d in 2:D)
          A_est <- A_est + (1 / D) * A_hat[[d]]
    }
    else if (strategy == "exp"){
      
      A <- error_mat
      beta <- median(apply(A, 2, function(x) sum(x ^ 2)))
      ## beta <- min(apply(A, 2, function(x) sum(x ^ 2))) / 10
      ## beta <- max(apply(A, 2, function(x) sum(x ^ 2))) * 2
      u_hat <- apply(A, 2, function(x) exp(sum(-x ^ 2 / beta)))
      u_hat <- u_hat / sum(u_hat)
      A_est <- u_hat[1] * A_hat[[1]]
      if(D != 1)
        for(d in 2:D)
          A_est <- A_est + u_hat[d] * A_hat[[d]]
    }
    else
      stop("strategy must be erm, average or exp!")

    As <- svd(A_est)
    K_hat <- As$u %*% diag(As$d / (1 - As$d)) %*% t(As$u)
    
    lambda0 <- tuning(y, K_hat, mode, lambda)
    K1 <- cbind(1, K_hat)
    K2 <- cbind(0, rbind(0, K_hat))
    
    theta <- ginv(lambda0 * K2 + t(K1) %*% K1) %*% t(K1) %*% y
    beta0 <- theta[1]
    alpha <- theta[-1]
    
    return(list(lam = lambda0, intercept = beta0, 
                alpha = alpha, K = K_hat, u_hat = u_hat))
  }


baseEstimate <- function(size, magn, Y, X1, X2, Kernlist, mode, lambda){
  
  A_hat <- list()
  error_mat <- matrix(0, nrow = size, ncol = magn)
  
  for (d in seq(magn)){
    Kern <- Kernlist[[d]]
    K1_m <- Kern(X1, X1)
    K2_m <- Kern(X2, X2)
    if(tr(K1_m) > 0 & tr(K2_m) > 0){
      K1_m <- K1_m / tr(K1_m)
      K2_m <- K2_m / tr(K2_m)
    }
    K <- K1_m + K2_m
    if (length(lambda) != 1){
      lambda0 <- tuning(Y, K, mode, lambda)
      K1 <- cbind(1, K)
      K2 <- cbind(0, rbind(0, K))
      theta <- ginv(lambda0 * K2 + t(K1) %*% K1) %*% t(K1) %*% Y
      beta0 <- theta[1]
      M <- K %*% ginv(K + lambda0 * diag(size))
      error_mat[, d] <- (diag(size) - M) %*% (Y - beta0) / diag(diag(size) - M)
      A_hat[[d]] <- M
    }
  }
  
  return(list(A_hat = A_hat, error_mat = error_mat))
}
