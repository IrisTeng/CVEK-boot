Ensemble <- 
  function(formula, label.names, Kernlist, data = NULL, beta = 1, 
           lambda = exp(seq(-5, 5, 1)), mode = 'loocv', strategy = 'erm'){
    
    y <- data[, as.character(attr(terms(formula), "variables"))[2]]
    
    # extract the information from the given formula to construct a true formula
    re <- GenericFormula(formula, label.names)
    generic.formula0 <- re[[1]]
    len <- re[[2]]
    
    # extract design matrix
    X <- model.matrix(generic.formula0, data)[, -1]
    n <- nrow(X)
    
    # estimation
    # D = 10
    X1 <- X[, c(1:length(label.names[[1]]))]
    X2 <- X[, c((length(label.names[[1]]) + 1):len)]
    
    D <- length(Kernlist)
    A.hat <- list()
    # lambda.mat <- matrix(0, nrow = D, ncol = n + 1)
    error.mat <- matrix(0, nrow = n, ncol = D)
    
    for (d in seq(D)){
      Kern <- Kernlist[[d]]
      K1.m <- Kern(X1, X1)
      K2.m <- Kern(X2, X2)
	  if(tr(K1.m) > 0 & tr(K2.m) > 0){
	    K1.m <- K1.m / tr(K1.m)
        K2.m <- K2.m / tr(K2.m)
	  }
      K <- K1.m + K2.m
      if (length(lambda) != 1){
        lambda0 <- Tuning(y, K, lambda, mode = mode)
        K1 <- cbind(1, K)
        K2 <- cbind(0, rbind(0, K))
        theta <- ginv(lambda0 * K2 + t(K1) %*% K1) %*% t(K1) %*% y
        # error.mat[, d] <- y - K1 %*% theta
        beta0 <- theta[1]
        M <- K %*% ginv(K + lambda0 * diag(n))
        error.mat[, d] <- (diag(n) - M) %*% (y - beta0) / diag(diag(n) - M)
        A.hat[[d]] <- M
      }
    }
    
    if (strategy == "erm"){
      
      A <- error.mat
      B <- rep(0, n)
      E <- rep(1, D)
      F <- 1
      G <- diag(D)
      H <- rep(0, D)
      u.hat <- lsei(A, B, E = E, F = F, G = G, H = H)$X
      A.est <- u.hat[1] * A.hat[[1]]
      if(D != 1)
        for(d in 2:D)
          A.est <- A.est + u.hat[d] * A.hat[[d]]
    }
    else if (strategy == "average"){
      
	  u.hat <- rep(1/D, D)
      A.est <- (1 / D) * A.hat[[1]]
      if(D != 1)
        for(d in 2:D)
          A.est <- A.est + (1 / D) * A.hat[[d]]
    }
    else if (strategy == "exp"){
      
      A <- error.mat
      u.hat <- apply(A, 2, function(x) exp(sum(-x ^ 2 / beta)))
      u.hat <- u.hat / sum(u.hat)
      A.est <- u.hat[1] * A.hat[[1]]
      if(D != 1)
        for(d in 2:D)
          A.est <- A.est + u.hat[d] * A.hat[[d]]
    }
    else
      stop("strategy must be loocv or average!")

    As <- svd(A.est)
    K.hat <- As$u %*% diag(As$d / (1 - As$d)) %*% t(As$u)
    
    lambda0 <- Tuning(y, K.hat, lambda, mode = mode)
    K1 <- cbind(1, K.hat)
    K2 <- cbind(0, rbind(0, K.hat))
    
    theta <- ginv(lambda0 * K2 + t(K1) %*% K1) %*% t(K1) %*% y
    beta0 <- theta[1]
    alpha <- theta[-1]
    
    return(list(sigma2.n = lambda0, intercept = beta0, 
                alpha = alpha, K = K.hat, u.hat = u.hat))
  }
