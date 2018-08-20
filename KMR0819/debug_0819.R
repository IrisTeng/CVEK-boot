library('devtools') 
setwd('/Users/dorabeedeng/Desktop')
setwd('./CVEK')
load_all()

label_names <- list(X1 = c("x1", "x2"), X2 = c("x3", "x4"))

set.seed(917)
data <- generate_data(n = 50, label_names, method = "rbf", 
                      int_effect = .3, l = 1, eps = .01)

kern_par <- data.frame(method = "linear", 
                       Sigma = 0, l = 1, p = 1)
kern_par$method <- as.character(kern_par$method)

kern <- generate_kernel("linear", p = 1)

formula <- Y ~ X1 + X2
fit <- define_model(formula, label_names, data, kern_par)
mode <- "loocv"
strategy <- "avg"
lambda <- exp(seq(-5, 5))

sol <- estimation(fit$Y, fit$X1, fit$X2, fit$kern_list, mode, strategy, lambda)
yhat_kern <- sol$intercept + sol$K %*% sol$alpha


set.seed(6370)
data_test <- generate_data(n = 40, label_names, method = "rbf", 
                           int_effect = .3, l = 1, eps = .01)
fit_test <- define_model(formula, label_names, data_test, kern_par)

formula_int <- Y ~ X1 * X2
test <- "boot"
B <- 100

## using CVEK
out_kern <- testing(formula_int, label_names, fit$Y, fit$X1, fit$X2, fit$kern_list, 
                    mode, strategy, beta = 1, test, lambda, B, data_test, fit_test)


## using lm.ridge
mod_train <- lm.ridge(Y ~ ., data, lambda = exp(seq(-5, 5)))
min_train <- which.min(mod_train$GCV)
coef_train <- coef(mod_train)[min_train, ]
X_train <- data[, -1]
X_train <- cbind(1, X_train)
X_train <- as.matrix(X_train)
y_train <- X_train %*% coef_train

yhat_comp <- cbind(yhat_kern, y_train, data$Y)
colnames(yhat_comp) <- c("CVEK", "lmridge", "true")
View(yhat_comp)


## look into function estimation
Y <- fit$Y 
X1 <- fit$X1 
X2 <- fit$X2
kern_list <- fit$kern_list

n <- length(Y)
kern_size <- length(kern_list)
out <- estimate_base(n, kern_size, Y, X1, X2, kern_list, mode, lambda)
A_hat <- out$A_hat
error_mat <- out$error_mat

out2 <- ensemble(n, kern_size, strategy, beta, error_mat, A_hat)
A_est <- out2$A_est
u_hat <- out2$u_hat

As <- svd(A_est)
K_hat <- As$u %*% diag(As$d / (1 - As$d)) %*% t(As$u)

####################################################################
## see if the same 
K_linear <- kern(X1, X1) / tr(kern(X1, X1)) + kern(X2, X2) / tr(kern(X2, X2))
K_linear / K_hat
## all equal to exp(5)

## check kernel implementation
K_imp <- X1 %*% t(X1) / tr(X1 %*% t(X1)) + X2 %*% t(X2) / tr(X2 %*% t(X2))
K_imp / K_hat
## all equal to exp(5), implementation is correct
####################################################################

## go on, check function tuning
lambda0 <- tuning(Y, K_hat, mode, lambda)
lambda0 == min_train
## not the same!!!  (>_<)

## just keep going for now!
K1 <- cbind(1, K_hat)
K2 <- cbind(0, rbind(0, K_hat))

theta <- ginv(lambda0 * K2 + t(K1) %*% K1) %*% t(K1) %*% Y
beta0 <- theta[1]
alpha <- theta[-1]

y_temp <- beta0 + K_hat %*% alpha
yhat_temp <- cbind(y_temp, yhat_comp)
colnames(yhat_temp)[1] <- "CVEK_temp"
View(yhat_temp)
## the first two columns are the same

## what if we condition on the same lambda=min_train?
lambda0 <- min_train
K1 <- cbind(1, K_hat)
K2 <- cbind(0, rbind(0, K_hat))

theta <- ginv(lambda0 * K2 + t(K1) %*% K1) %*% t(K1) %*% Y
beta0 <- theta[1]
alpha <- theta[-1]

y_con <- beta0 + K_hat %*% alpha
yhat_con <- cbind(y_con, yhat_temp)
colnames(yhat_con)[1] <- "CVEK_con"
View(yhat_con)
## ao! estimation is dangerous! (>=<)


