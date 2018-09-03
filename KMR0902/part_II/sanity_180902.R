library("devtools") 
install_github("IrisTeng/CVEK")
library(CVEK)

label_names <- list(X1 = c("x1", "x2"), X2 = c("x3", "x4"))
formula <- Y ~ X1 + X2
formula_int <- Y ~ X1 * X2
method <- "polynomial"
int_effect <- 0
eps <- .01
beta <- "min"
test <- "boot"
lambda <- exp(seq(-5, 5, .2))
B <- 100

for (p in 1:3) {
  l <- p / 2
  kern <- data.frame(method = method, Sigma = 0, l = l, p = p)
  kern$method <- as.character(kern$method)
  data_all <- generate_data(200, label_names, method = method,
                            int_effect = int_effect, l = l, 
                            p = p, eps = eps)
  data_train <- data_all[1:100, ]
  data_test <- data_all[101:200, ]
  fit <- define_model(formula, label_names, data_train, kern)
  fit_test <- define_model(formula, label_names, data_test, kern)
  for (mode in c("AICc", "GCVc", "gmpml", "loocv")) {
    for (strategy in c("avg", "exp", "erm")) {
      res_tst <- testing(formula_int, label_names, fit$Y, fit$X1, fit$X2, 
                         fit$kern_list, mode, strategy, beta, test, lambda, 
                         B, data_test, fit_test)
      msg <- paste0("true=", method, ", p=", p, ", mode=", 
                    mode, ", strategy=", strategy)
      output <- paste0("selected lambda=", round(log(res_tst$lam), 2), 
                       ", train_RMSE=", round(res_tst$train_RMSE, 4), 
                       ", test_RMSE=", round(res_tst$test_RMSE, 4))
      print(msg)
      print(output)
    }
  }
}


