tune_ridge <- function(data, lambda) {
  
  n <- nrow(data)
  Y0 <- data[, 1]
  beta0 <- mean(Y0)
  Y <- Y0 - beta0
  X <- as.matrix(data[, -1])
  p <- ncol(X)
  Xm <- colMeans(X)
  X <- X - rep(Xm, rep(n, p))
  Xscale <- drop(rep(1 / n, n) %*% X ^ 2) ^ .5
  X <- X / rep(Xscale, rep(n, p))
  CV <- sapply(lambda, function(k) {
    
    A <- X %*% ginv(t(X) %*% X + k * diag(p)) %*% t(X)
    sum(((diag(n) - A) %*% Y / diag(diag(n) - A)) ^ 2)
  })
  
  lambda0 <- lambda[which(CV == min(CV))]
  Y_hat <- X %*% ginv(t(X) %*% X + lambda0 * diag(p)) %*% t(X) %*% Y + beta0
  
  c(beta0, ginv(t(X) %*% X + lambda0 * diag(p)) %*% t(X) %*% Y)
}


cet_X <- function(data) {
  
  n <- nrow(data)
  X <- as.matrix(data[, -1])
  p <- ncol(X)
  Xm <- colMeans(X)
  X <- X - rep(Xm, rep(n, p))
  Xscale <- drop(rep(1 / n, n) %*% X ^ 2) ^ .5
  X <- X / rep(Xscale, rep(n, p))
  cbind(1, X)
}

beta_man <- tune_ridge(data, lambda)
yman_hat <- cet_X(data) %*% beta_man
yman_test <- cet_X(data_test) %*% beta_man

ss <- cbind(data$Y, out$y_hat, yman_hat, data_test$Y, 
            out$y_test, yman_test)
colnames(ss) <- c("train_true", "train_kern", "train_rid", 
                  "test_true", "test_kern", "test_rid")




