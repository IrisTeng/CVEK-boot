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


# testCV <-
#   function(y, X, Kern, lambda = exp(seq(-5, 5, 1))){
#     n <- nrow(X)
#     cv_idx <- sample(1:n, round(n/2))
#     group_id <- list(1:2, 3:4)
# 
#     y_tr <- y[-cv_idx]; y_cv <- y[cv_idx]
#     X_tr <- X[-cv_idx, ]; X_cv <- X[cv_idx, ]
#     K_tr <-
#       (Kern(X_tr[, group_id[[1]]], X_tr[, group_id[[1]]]) +
#          Kern(X_tr[, group_id[[2]]], X_tr[, group_id[[2]]])
#       )/length(group_id)
#     K_cv <-
#       (Kern(X_cv[, group_id[[1]]], X_tr[, group_id[[1]]]) +
#          Kern(X_cv[, group_id[[2]]], X_tr[, group_id[[2]]])
#       )/length(group_id)
# 
# 
#     # estimation
#     CV <- lapply(lambda, function(k){
#       # train
#       # K1 <- cbind(1, K_tr)
#       # K2 <- cbind(0, rbind(0, K_tr))
#       #
#       # theta <- ginv(k * K2 + t(K1) %*% K1) %*% t(K1) %*% y_tr
#       # beta_tr <- theta[1]
#       # alpha_tr <- theta[-1]
#       alpha_tr <- ginv(k * diag(nrow(K_tr)) + K_tr) %*% y_tr
# 
#       # train and test
#       y_pred_tr <- K_tr %*% alpha_tr
#       y_pred_cv <- K_cv %*% alpha_tr
# 
#       c(sum((y_tr - y_pred_tr)^2)/length(y_tr),
#         sum((y_cv - y_pred_cv)^2)/length(y_cv))
#     })
# 
#     CV <- do.call("rbind", CV)
# 
#     plot(rep(order(lambda), 2), log(CV), type = "n")
#     lines(order(lambda), log(CV[, 1]), col = 1)
#     lines(order(lambda), log(CV[, 2]), col = 2)
# 
# 
#     lambda0 <- lambda[which(CV[, 2] == min(CV[, 2]))]
#     return(lambda0)
#   }
