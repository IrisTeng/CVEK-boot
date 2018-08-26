
produce_design_matrix <- function(X1, X2, kern){
  U_1 <- svd(kern(X1, X1))$u
  d_1 <- svd(kern(X1, X1))$d
  X1_pseudo <- U_1 %*% sqrt(diag(d_1))[, which(d_1 > 1e-12)]
  
  U_2 <- svd(kern(X2, X2))$u
  d_2 <- svd(kern(X2, X2))$d
  X2_pseudo <- U_2 %*% sqrt(diag(d_2))[, which(d_2 > 1e-12)]
  
  cbind(1, X1_pseudo, X2_pseudo)
}

label_names <- list(X1 = c("x1", "x2"), X2 = c("x3", "x4"))
mode <- "loocv"
strategy <- "avg"
lambda <- .5
n <- 200

####################################################################
## linear
set.seed(6370)
data_all <- generate_data(n = 200, label_names, method = "linear", 
                          int_effect = 0, l = 1, p = 1, eps = .01)
kern <- generate_kernel("linear", p = 1)

data <- data_all[1:100,]
data_test <- data_all[101:200,]

Y <- data$Y
Y_new <- data_test$Y
X1 <- as.matrix(data[, 2:3])
X2 <- as.matrix(data[, 4:5])
X <- produce_design_matrix(X1, X2, kern)

X1_new <- as.matrix(data_test[, 2:3])
X2_new <- as.matrix(data_test[, 4:5])
X_new <- produce_design_matrix(X1_new, X2_new, kern)

p <- ncol(X) - 1
# ridge regression should be like this:
I <- cbind(0, rbind(0, diag(p)))
beta_ridge <- ginv(t(X) %*% X + lambda * I) %*% t(X) %*% Y
pred_train_ridge <- X %*% beta_ridge
pred_test_ridge <- X_new %*% beta_ridge
train_RMSE0 <- sqrt(sum((pred_train_ridge - Y) ^ 2) / n)
test_RMSE0 <- sqrt(sum((pred_test_ridge - Y_new) ^ 2) / n)
train_RMSE0 / test_RMSE0

## ridge again
X <- cbind(1, X1, X2)
X_new <- cbind(1, X1_new, X2_new)
p <- ncol(X) - 1
# ridge regression should be like this:
I <- cbind(0, rbind(0, diag(p)))
beta_ridge <- ginv(t(X) %*% X + lambda * I) %*% t(X) %*% Y
pred_train_ridge <- X %*% beta_ridge
pred_test_ridge <- X_new %*% beta_ridge
train_RMSE0 <- sqrt(sum((pred_train_ridge - Y) ^ 2) / n)
test_RMSE0 <- sqrt(sum((pred_test_ridge - Y_new) ^ 2) / n)
train_RMSE0 / test_RMSE0


# kernel regression should be like this:
K <- kern(X1, X1) + kern(X2, X2)
K_test <- kern(X1_new, X1) + kern(X2_new, X2)
# K <- X[, -1] %*% t(X[, -1])
# K_test <- X_new[, -1] %*% t(X[, -1])
K1_test <- cbind(1, K_test)
K1 <- cbind(1, K)
K2 <- cbind(0, rbind(0, K))
theta <- ginv(lambda * K2 + t(K1) %*% K1) %*% t(K1) %*% Y
beta0_lkr <- theta[1]
alpha_lkr <- theta[-1]
pred_train_lkr <- beta0_lkr + K %*% alpha_lkr
pred_test_lkr <- beta0_lkr + K_test %*% alpha_lkr
plot(pred_train_lkr, pred_train_ridge)
pred_test_lkr <- beta0_lkr + K_test %*% alpha_lkr
train_RMSE <- sqrt(sum((pred_train_lkr - Y) ^ 2) / n)
test_RMSE <- sqrt(sum((pred_test_lkr - Y_new) ^ 2) / n)
train_RMSE / test_RMSE
####################################################################


####################################################################
## quadratic
set.seed(917)
data_all <- generate_data(n = 200, label_names, method = "polynomial", 
                          int_effect = 0, l = 1, p = 2, eps = .01)
kern <- generate_kernel("polynomial", p = 2)

data <- data_all[1:100,]
data_test <- data_all[101:200,]

Y <- data$Y
Y_new <- data_test$Y
X1 <- as.matrix(data[, 2:3])
X2 <- as.matrix(data[, 4:5])
X <- produce_design_matrix(X1, X2, kern)

X1_new <- as.matrix(data_test[, 2:3])
X2_new <- as.matrix(data_test[, 4:5])
X_new <- produce_design_matrix(X1_new, X2_new, kern)

p <- ncol(X) - 1
# ridge regression should be like this:
I <- cbind(0, rbind(0, diag(p)))
beta_ridge <- ginv(t(X) %*% X + lambda * I) %*% t(X) %*% Y
pred_train_ridge <- X %*% beta_ridge
pred_test_ridge <- X_new %*% beta_ridge
train_RMSE0 <- sqrt(sum((pred_train_ridge - Y) ^ 2) / n)
test_RMSE0 <- sqrt(sum((pred_test_ridge - Y_new) ^ 2) / n)
train_RMSE0 / test_RMSE0

## ridge again
X <- cbind(1, X1, X1 ^ 2, X2, X2 ^ 2)
X_new <- cbind(1, X1_new, X1_new ^ 2, X2_new, X2_new ^ 2)
p <- ncol(X) - 1
# ridge regression should be like this:
I <- cbind(0, rbind(0, diag(p)))
beta_ridge <- ginv(t(X) %*% X + lambda * I) %*% t(X) %*% Y
pred_train_ridge <- X %*% beta_ridge
pred_test_ridge <- X_new %*% beta_ridge
train_RMSE0 <- sqrt(sum((pred_train_ridge - Y) ^ 2) / n)
test_RMSE0 <- sqrt(sum((pred_test_ridge - Y_new) ^ 2) / n)
train_RMSE0 / test_RMSE0


# kernel regression should be like this:
K <- kern(X1, X1) + kern(X2, X2)
K_test <- kern(X1_new, X1) + kern(X2_new, X2)
# K <- X[, -1] %*% t(X[, -1])
# K_test <- X_new[, -1] %*% t(X[, -1])
K1_test <- cbind(1, K_test)
K1 <- cbind(1, K)
K2 <- cbind(0, rbind(0, K))
theta <- ginv(lambda * K2 + t(K1) %*% K1) %*% t(K1) %*% Y
beta0_lkr <- theta[1]
alpha_lkr <- theta[-1]
pred_train_lkr <- beta0_lkr + K %*% alpha_lkr
pred_test_lkr <- beta0_lkr + K_test %*% alpha_lkr
plot(pred_train_lkr, pred_train_ridge)
pred_test_lkr <- beta0_lkr + K_test %*% alpha_lkr
train_RMSE <- sqrt(sum((pred_train_lkr - Y) ^ 2) / n)
test_RMSE <- sqrt(sum((pred_test_lkr - Y_new) ^ 2) / n)
train_RMSE / test_RMSE
####################################################################


####################################################################
## cubic
set.seed(0212)
data_all <- generate_data(n = 200, label_names, method = "polynomial", 
                          int_effect = 0, l = 1, p = 3, eps = .01)
kern <- generate_kernel("polynomial", p = 3)

data <- data_all[1:100,]
data_test <- data_all[101:200,]

Y <- data$Y
Y_new <- data_test$Y
X1 <- as.matrix(data[, 2:3])
X2 <- as.matrix(data[, 4:5])
X <- produce_design_matrix(X1, X2, kern)

X1_new <- as.matrix(data_test[, 2:3])
X2_new <- as.matrix(data_test[, 4:5])
X_new <- produce_design_matrix(X1_new, X2_new, kern)

p <- ncol(X) - 1
# ridge regression should be like this:
I <- cbind(0, rbind(0, diag(p)))
beta_ridge <- ginv(t(X) %*% X + lambda * I) %*% t(X) %*% Y
pred_train_ridge <- X %*% beta_ridge
pred_test_ridge <- X_new[, -22] %*% beta_ridge
train_RMSE0 <- sqrt(sum((pred_train_ridge - Y) ^ 2) / n)
test_RMSE0 <- sqrt(sum((pred_test_ridge - Y_new) ^ 2) / n)
train_RMSE0 / test_RMSE0

## ridge again
X <- cbind(1, X1, X1 ^ 2, X1 ^ 3, X2, X2 ^ 2, X2 ^ 3)
X_new <- cbind(1, X1_new, X1_new ^ 2, X1_new ^ 3, X2_new, X2_new ^ 2,  X2_new ^ 3)
p <- ncol(X) - 1
# ridge regression should be like this:
I <- cbind(0, rbind(0, diag(p)))
beta_ridge <- ginv(t(X) %*% X + lambda * I) %*% t(X) %*% Y
pred_train_ridge <- X %*% beta_ridge
pred_test_ridge <- X_new %*% beta_ridge
train_RMSE0 <- sqrt(sum((pred_train_ridge - Y) ^ 2) / n)
test_RMSE0 <- sqrt(sum((pred_test_ridge - Y_new) ^ 2) / n)
train_RMSE0 / test_RMSE0


# kernel regression should be like this:
K <- kern(X1, X1) + kern(X2, X2)
K_test <- kern(X1_new, X1) + kern(X2_new, X2)
# K <- X[, -1] %*% t(X[, -1])
# K_test <- X_new[, -1] %*% t(X[, -1])
K1_test <- cbind(1, K_test)
K1 <- cbind(1, K)
K2 <- cbind(0, rbind(0, K))
theta <- ginv(lambda * K2 + t(K1) %*% K1) %*% t(K1) %*% Y
beta0_lkr <- theta[1]
alpha_lkr <- theta[-1]
pred_train_lkr <- beta0_lkr + K %*% alpha_lkr
pred_test_lkr <- beta0_lkr + K_test %*% alpha_lkr
plot(pred_train_lkr, pred_train_ridge)
pred_test_lkr <- beta0_lkr + K_test %*% alpha_lkr
train_RMSE <- sqrt(sum((pred_train_lkr - Y) ^ 2) / n)
test_RMSE <- sqrt(sum((pred_test_lkr - Y_new) ^ 2) / n)
train_RMSE / test_RMSE
####################################################################


####################################################################
## rbf
set.seed(4486)
data_all <- generate_data(n = 200, label_names, method = "rbf", 
                          int_effect = 0, l = 1, p = 2, eps = .01)
kern <- generate_kernel("rbf", l = 1)

data <- data_all[1:100,]
data_test <- data_all[101:200,]

Y <- data$Y
Y_new <- data_test$Y
X1 <- as.matrix(data[, 2:3])
X2 <- as.matrix(data[, 4:5])
X <- produce_design_matrix(X1, X2, kern)[, -c(186:189)]

X1_new <- as.matrix(data_test[, 2:3])
X2_new <- as.matrix(data_test[, 4:5])
X_new <- produce_design_matrix(X1_new, X2_new, kern)

p <- ncol(X) - 1
# ridge regression should be like this:
I <- cbind(0, rbind(0, diag(p)))
beta_ridge <- ginv(t(X) %*% X + lambda * I) %*% t(X) %*% Y
pred_train_ridge <- X %*% beta_ridge
pred_test_ridge <- X_new %*% beta_ridge
train_RMSE0 <- sqrt(sum((pred_train_ridge - Y) ^ 2) / n)
test_RMSE0 <- sqrt(sum((pred_test_ridge - Y_new) ^ 2) / n)
train_RMSE0 / test_RMSE0

## ridge again
X <- cbind(1, X1, X1 ^ 2, X1 ^ 3, X2, X2 ^ 2, X2 ^ 3)
X_new <- cbind(1, X1_new, X1_new ^ 2, X1_new ^ 3, X2_new, X2_new ^ 2,  X2_new ^ 3)
p <- ncol(X) - 1
# ridge regression should be like this:
I <- cbind(0, rbind(0, diag(p)))
beta_ridge <- ginv(t(X) %*% X + lambda * I) %*% t(X) %*% Y
pred_train_ridge <- X %*% beta_ridge
pred_test_ridge <- X_new %*% beta_ridge
train_RMSE0 <- sqrt(sum((pred_train_ridge - Y) ^ 2) / n)
test_RMSE0 <- sqrt(sum((pred_test_ridge - Y_new) ^ 2) / n)
train_RMSE0 / test_RMSE0


# kernel regression should be like this:
K <- kern(X1, X1) + kern(X2, X2)
K_test <- kern(X1_new, X1) + kern(X2_new, X2)
# K <- X[, -1] %*% t(X[, -1])
# K_test <- X_new[, -1] %*% t(X[, -1])
K1_test <- cbind(1, K_test)
K1 <- cbind(1, K)
K2 <- cbind(0, rbind(0, K))
theta <- ginv(lambda * K2 + t(K1) %*% K1) %*% t(K1) %*% Y
beta0_lkr <- theta[1]
alpha_lkr <- theta[-1]
pred_train_lkr <- beta0_lkr + K %*% alpha_lkr
pred_test_lkr <- beta0_lkr + K_test %*% alpha_lkr
plot(pred_train_lkr, pred_train_ridge)
pred_test_lkr <- beta0_lkr + K_test %*% alpha_lkr
train_RMSE <- sqrt(sum((pred_train_lkr - Y) ^ 2) / n)
test_RMSE <- sqrt(sum((pred_test_lkr - Y_new) ^ 2) / n)
train_RMSE / test_RMSE
####################################################################
