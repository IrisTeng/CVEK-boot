gen_data <- 
  function(n = 100, p = 3, int_str = 0){
    X <- matrix(rnorm(n = 2*n*p), nrow = n)
    # produce interaction term
    X1 <- X[, 1:(p)]
    X2 <- X[, (1 + p):(2 * p)]
    
    X12 <- NULL
    for (i in 1:p){
      X12 <- cbind(X12, X1[, i] * X2)
    }
    
    #
    X <- cbind(1, X)
    P_X <- X %*% MASS::ginv(t(X) %*% X) %*% t(X)
    
    beta <- rnorm(n = 2*p + 1)
    beta_12 <- rnorm(n = p^2)
    
    h_0 <- X %*% beta
    h_12 <- X12 %*% beta_12
    h_0 <- h_0/sqrt(sum(h_0^2))
    h_12 <- h_12/sqrt(sum(h_12^2))
    
    eps <- rnorm(n = n, sd = 0.01)
    
    y <- h_0 + int_str * (diag(n) - P_X) %*% h_12 + eps
    
    list(y = y, X = X, beta = beta)
  }

test_stat <- 
  function(y, X, beta_0_est){
    X0 <- X[, -1]
    p <- ncol(X0)/2
    X1 <- X0[, 1:(p)]
    X2 <- X0[, (1 + p):(2 * p)]
    
    X12 <- NULL
    for (i in 1:p){
      X12 <- cbind(X12, X1[, i] * X2)
    }
    
    V_inv <- MASS::ginv(X %*% t(X))
    K12 <- X12 %*% t(X12)
    y_center <- matrix(y - beta_0_est, ncol = 1)
    
    t(y_center) %*% V_inv %*% K12 %*% V_inv %*% y_center
  }

#### test run ####

simu <- function(index, n = 100, int_str = 0, B = 200){
  data_list <- gen_data(n = n, int_str = int_str)
  
  # extract data and estimate
  X <- data_list$X
  y <- data_list$y
  lm <- lm(y ~ X)
  beta_0_est <- coef(lm)[1]
  y_est <- predict(lm)
  eps_est <- summary(lm(y ~ X))$sigma
  
  # compute statistic
  #index_boot <- sample(1:n, n, replace = TRUE)
  #X_boot <- X[index_boot, ]
  #y_boot <- y[index_boot] + rnorm(n, sd = eps_est)
  
  boot_stat <- sapply(1:B, function(k){
    y_boot <- y_est + rnorm(n, sd = eps_est)
    test_stat(y_boot, X, beta_0_est)
  })
  # y_boot <- y_est + rnorm(n, sd = eps_est)
    
  true_stat <- test_stat(y, X, beta_0_est)
  # boot_stat <- test_stat(y_boot, X, beta_0_est)
  
  pvalue <- sum(as.numeric(true_stat) <= boot_stat) / B
  
  return(list(pvalue = pvalue, boot.stat = boot_stat, true.stat = true_stat))
  # return(c(true_stat, boot_stat, sum((y - y_est)^2)))
}

output <- t(sapply(1, simu, n=100, int_str = 0))
# output2 <- t(sapply(1, simu, n=100, int_str = 0))

plot(density(output[, 1]))
lines(density(output[, 2]), col = 2)

