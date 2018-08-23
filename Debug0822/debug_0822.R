####################################################################
## linear
n <- 100
p <- 5
set.seed(0212)
X <- cbind(1, matrix(rnorm(n * p), nrow = n))
X_new <- cbind(1, matrix(rnorm(n * p), nrow = n))
beta <- rnorm(p + 1)
Y <- X %*% beta + rnorm(sd = .1, n = n)
Y_new <- X_new %*% beta + rnorm(sd = .1, n = n)
lambda <- .5

data <- cbind(Y, X[, -1])
X1 <- data[, 2:3]
X2 <- data[, 4:6]
data <- as.data.frame(data)
colnames(data) <- c("Y", paste0("x", 1:5))

# 
# mod <- lm.ridge(Y ~ ., data, lambda)
# min <- which.min(mod$GCV)
# coef_tr <- coef(mod)[min, ]

# Y <- data$Y
# X <- as.matrix(cbind(1, data[, -1]))
# p <- ncol(X) - 1
# ridge regression should be like this:
I <- cbind(0, rbind(0, diag(p)))
beta_ridge <- ginv(t(X) %*% X + lambda * I) %*% t(X) %*% Y
pred_train_ridge <- X %*% beta_ridge
pred_test_ridge <- X_new %*% beta_ridge

train_RMSE0 <- sqrt(sum((pred_train_ridge - Y) ^ 2) / n)
test_RMSE0 <- sqrt(sum((pred_test_ridge - Y_new) ^ 2) / n)
train_RMSE0 / test_RMSE0

# linear kernel regression should be like this:
K <- X[, -1] %*% t(X[, -1])
K_test <- X_new[, -1] %*% t(X[, -1])
K1_test <- cbind(1, K_test)
K1 <- cbind(1, K)
K2 <- cbind(0, rbind(0, K))

theta <- ginv(lambda * K2 + t(K1) %*% K1) %*% t(K1) %*% Y
beta0_lkr <- theta[1]
alpha_lkr <- theta[-1]
pred_train_lkr <- beta0_lkr + K %*% alpha_lkr
pred_test_lkr <- beta0_lkr + K_test %*% alpha_lkr

train_RMSE <- sqrt(sum((pred_train_lkr - Y) ^ 2) / n)
test_RMSE <- sqrt(sum((pred_test_lkr - Y_new) ^ 2) / n)
train_RMSE / test_RMSE


A <- K1 %*% ginv(t(K1) %*% K1 + lambda * K2) %*% t(K1)
A_test <- K1_test %*% ginv(t(K1) %*% K1 + lambda * K2) %*% t(K1)
pred_train_prj <- A %*% Y
pred_test_prj <- A_test %*% Y

train_RMSE <- sqrt(sum((pred_train_prj - Y) ^ 2) / n)
test_RMSE <- sqrt(sum((pred_test_prj - Y_new) ^ 2) / n)
train_RMSE / test_RMSE
####################################################################


####################################################################
## quadratic
label_names <- list(X1 = c("x1", "x2"), X2 = c("x3", "x4"))
mode <- "loocv"
strategy <- "avg"
lambda <- .5
kern <- generate_kernel("polynomial", p = 2)

set.seed(917)
data <- generate_data(n = 100, label_names, method = "polynomial", 
                      int_effect = 0, l = 1, p = 2, eps = .01)

set.seed(6370)
data_test <- generate_data(n = 100, label_names, method = "polynomial", 
                           int_effect = 0, l = 1, p = 2, eps = .01)

Y <- data$Y
Y_new <- data_test$Y
X1 <- data[, 2:3]
X2 <- data[, 4:5]
X <- cbind(1, svd(kern(X1, X1))$u %*% sqrt(svd(kern(X1, X1))$d), 
           svd(kern(X2, X2))$u %*% sqrt(svd(kern(X2, X2))$d))

X1_new <- data_test[, 2:3]
X2_new <- data_test[, 4:5]
X_new <- cbind(1, svd(kern(X1_new, X1_new))$u %*% sqrt(svd(kern(X1_new, X1_new))$d), 
               svd(kern(X2_new, X2_new))$u %*% sqrt(svd(kern(X2_new, X2_new))$d))

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
# K <- kern(X1, X1) + kern(X2, X2)
# K_test <- kern(X1_new, X1) + kern(X2_new, X2)
K <- X[, -1] %*% t(X[, -1])
K_test <- X_new[, -1] %*% t(X[, -1])
K1_test <- cbind(1, K_test)
K1 <- cbind(1, K)
K2 <- cbind(0, rbind(0, K))

theta <- ginv(lambda * K2 + t(K1) %*% K1) %*% t(K1) %*% Y
beta0_lkr <- theta[1]
alpha_lkr <- theta[-1]
pred_train_lkr <- beta0_lkr + K %*% alpha_lkr
pred_test_lkr <- beta0_lkr + K_test %*% alpha_lkr

train_RMSE <- sqrt(sum((pred_train_lkr - Y) ^ 2) / n)
test_RMSE <- sqrt(sum((pred_test_lkr - Y_new) ^ 2) / n)
train_RMSE / test_RMSE


A <- K1 %*% ginv(t(K1) %*% K1 + lambda * K2) %*% t(K1)
A_test <- K1_test %*% ginv(t(K1) %*% K1 + lambda * K2) %*% t(K1)
pred_train_prj <- A %*% Y
pred_test_prj <- A_test %*% Y

train_RMSE <- sqrt(sum((pred_train_prj - Y) ^ 2) / n)
test_RMSE <- sqrt(sum((pred_test_prj - Y_new) ^ 2) / n)
train_RMSE / test_RMSE

####################################################################


####################################################################
## cubic
label_names <- list(X1 = c("x1", "x2"), X2 = c("x3", "x4"))
mode <- "loocv"
strategy <- "avg"
lambda <- .5
kern <- generate_kernel("polynomial", p = 3)

set.seed(917)
data <- generate_data(n = 100, label_names, method = "polynomial", 
                      int_effect = 0, l = 1, p = 3, eps = .01)

set.seed(6370)
data_test <- generate_data(n = 100, label_names, method = "polynomial", 
                           int_effect = 0, l = 1, p = 3, eps = .01)

Y <- data$Y
Y_new <- data_test$Y
X1 <- data[, 2:3]
X2 <- data[, 4:5]
X <- cbind(1, svd(kern(X1, X1))$u %*% sqrt(svd(kern(X1, X1))$d), 
           svd(kern(X2, X2))$u %*% sqrt(svd(kern(X2, X2))$d))

X1_new <- data_test[, 2:3]
X2_new <- data_test[, 4:5]
X_new <- cbind(1, svd(kern(X1_new, X1_new))$u %*% sqrt(svd(kern(X1_new, X1_new))$d), 
               svd(kern(X2_new, X2_new))$u %*% sqrt(svd(kern(X2_new, X2_new))$d))

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
# K <- kern(X1, X1) + kern(X2, X2)
# K_test <- kern(X1_new, X1) + kern(X2_new, X2)
K <- X[, -1] %*% t(X[, -1])
K_test <- X_new[, -1] %*% t(X[, -1])
K1_test <- cbind(1, K_test)
K1 <- cbind(1, K)
K2 <- cbind(0, rbind(0, K))

theta <- ginv(lambda * K2 + t(K1) %*% K1) %*% t(K1) %*% Y
beta0_lkr <- theta[1]
alpha_lkr <- theta[-1]
pred_train_lkr <- beta0_lkr + K %*% alpha_lkr
pred_test_lkr <- beta0_lkr + K_test %*% alpha_lkr

train_RMSE <- sqrt(sum((pred_train_lkr - Y) ^ 2) / n)
test_RMSE <- sqrt(sum((pred_test_lkr - Y_new) ^ 2) / n)
train_RMSE / test_RMSE


A <- K1 %*% ginv(t(K1) %*% K1 + lambda * K2) %*% t(K1)
A_test <- K1_test %*% ginv(t(K1) %*% K1 + lambda * K2) %*% t(K1)
pred_train_prj <- A %*% Y
pred_test_prj <- A_test %*% Y

train_RMSE <- sqrt(sum((pred_train_prj - Y) ^ 2) / n)
test_RMSE <- sqrt(sum((pred_test_prj - Y_new) ^ 2) / n)
train_RMSE / test_RMSE

####################################################################



####################################################################
## rbf
label_names <- list(X1 = c("x1", "x2"), X2 = c("x3", "x4"))
mode <- "loocv"
strategy <- "avg"
lambda <- .5
kern <- generate_kernel("rbf", l = 1)

set.seed(917)
data <- generate_data(n = 100, label_names, method = "rbf", 
                      int_effect = 0, l = 1, p = 2, eps = .01)

set.seed(6370)
data_test <- generate_data(n = 100, label_names, method = "rbf", 
                           int_effect = 0, l = 1, p = 2, eps = .01)

Y <- data$Y
Y_new <- data_test$Y
X1 <- data[, 2:3]
X2 <- data[, 4:5]
X <- cbind(1, svd(kern(X1, X1))$u %*% sqrt(svd(kern(X1, X1))$d), 
           svd(kern(X2, X2))$u %*% sqrt(svd(kern(X2, X2))$d))

X1_new <- data_test[, 2:3]
X2_new <- data_test[, 4:5]
X_new <- cbind(1, svd(kern(X1_new, X1_new))$u %*% sqrt(svd(kern(X1_new, X1_new))$d), 
               svd(kern(X2_new, X2_new))$u %*% sqrt(svd(kern(X2_new, X2_new))$d))

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
# K <- kern(X1, X1) + kern(X2, X2)
# K_test <- kern(X1_new, X1) + kern(X2_new, X2)
K <- X[, -1] %*% t(X[, -1])
K_test <- X_new[, -1] %*% t(X[, -1])
K1_test <- cbind(1, K_test)
K1 <- cbind(1, K)
K2 <- cbind(0, rbind(0, K))

theta <- ginv(lambda * K2 + t(K1) %*% K1) %*% t(K1) %*% Y
beta0_lkr <- theta[1]
alpha_lkr <- theta[-1]
pred_train_lkr <- beta0_lkr + K %*% alpha_lkr
pred_test_lkr <- beta0_lkr + K_test %*% alpha_lkr

train_RMSE <- sqrt(sum((pred_train_lkr - Y) ^ 2) / n)
test_RMSE <- sqrt(sum((pred_test_lkr - Y_new) ^ 2) / n)
train_RMSE / test_RMSE


A <- K1 %*% ginv(t(K1) %*% K1 + lambda * K2) %*% t(K1)
A_test <- K1_test %*% ginv(t(K1) %*% K1 + lambda * K2) %*% t(K1)
pred_train_prj <- A %*% Y
pred_test_prj <- A_test %*% Y

train_RMSE <- sqrt(sum((pred_train_prj - Y) ^ 2) / n)
test_RMSE <- sqrt(sum((pred_test_prj - Y_new) ^ 2) / n)
train_RMSE / test_RMSE

####################################################################
