#' Conducting Score Tests for Interaction
#'
#' Conduct score tests comparing a fitted model and a more general alternative
#' model.
#'
#' There are two tests available here:
#'
#' \bold{Asymptotic Test}
#'
#' This is based on the classical variance component test to construct a
#' testing procedure for the hypothesis about Gaussian process function.
#'
#' \bold{Bootstrap Test}
#'
#' When it comes to small sample size, we can use bootstrap test instead, which
#' can give valid tests with moderate sample sizes and requires similar
#' computational effort to a permutation test.
#'
#' @param formula_int (formula) A symbolic description of the model with 
#' interaction.
#' @param label_names (list) A character string indicating all the interior 
#' variables included in each predictor.
#' @param Y (vector of length n) Reponses of the dataframe.
#' @param X1 (dataframe, n*p1) The first type of factor in the dataframe (could 
#' contains several subfactors).
#' @param X2 (dataframe, n*p2) The second type of factor in the dataframe (could 
#' contains several subfactors).
#' @param kern_list (list of length K) A list of kernel functions given by user.
#' @param mode (character) A character string indicating which tuning 
#' parameter criteria is to be used.
#' @param strategy (character) A character string indicating which ensemble 
#' strategy is to be used.
#' @param beta (numeric/character) A numeric value specifying the parameter 
#' when strategy = "exp" \code{\link{ensemble_exp}}.
#' @param test (character) A character string indicating which test is 
#' to be used.
#' @param lambda (numeric) A numeric string specifying the range of 
#' noise to be chosen. The lower limit of lambda must be above 0.
#' @param B (integer) A numeric value indicating times of resampling w
#' hen test = "boot".
#' @return \item{pvalue}{(numeric) p-value of the test.}
#' \item{u_weight}{(numeric) p-value of the test.}
#' \item{lam}{(numeric) p-value of the test.}
#' \item{train_RMSE}{(numeric) p-value of the test.}
#' \item{test_RMSE}{(numeric) p-value of the test.}
#' \item{K_tr}{(numeric) p-value of the test.}
#' \item{A_tr}{(numeric) p-value of the test.}
#' \item{V0_inv_tr}{(numeric) p-value of the test.}
#' \item{K_eig}{(numeric) p-value of the test.}
#' @author Wenying Deng
#' @seealso method: \code{\link{generate_kernel}}
#'
#' mode: \code{\link{tuning}}
#'
#' strategy: \code{\link{ensemble}}
#' @references Xihong Lin. Variance component testing in generalised linear
#' models with random effects. June 1997.
#'
#' Arnab Maity and Xihong Lin. Powerful tests for detecting a gene effect in
#' the presence of possible gene-gene interactions using garrote kernel
#' machines. December 2011.
#'
#' Petra Bu ̊zˇkova ́, Thomas Lumley, and Kenneth Rice. Permutation and
#' parametric bootstrap tests for gene-gene and gene-environment interactions.
#' January 2011.
#' @examples
#'
#'
#' testing(formula_int = Y ~ X1 * X2,
#' label_names = list(X1 = c("x1", "x2"), X2 = c("x3", "x4")),
#' Y, X1, X2, kern_list, mode = "loocv", strategy = "erm",
#' beta = 1, test = "boot", lambda = exp(seq(-5, 5)), B = 100)
#'
#'
#' @export testing

testing <- function(formula_int, label_names, Y, X1, X2, kern_list,
                    mode = "loocv", strategy = "erm", beta = 1,
                    test = "boot", lambda_list = exp(seq(-5, 5)), 
                    B = 100, data_test = NULL, fit_test = NULL) {
  
  re <- generate_formula(formula_int, label_names)
  generic_formula0 <- re$generic_formula
  len <- re$length_main
  data <- as.data.frame(cbind(Y, X1, X2))
  colnames(data) <- c("Y", label_names[[1]], label_names[[2]])
  X <- model.matrix(generic_formula0, data)[, -1]
  X12 <- X[, c((len + 1):dim(X)[2])]
  n <- length(Y)
  
  result <- estimation(Y, X1, X2, kern_list, mode, strategy, beta, lambda_list)
  lambda <- result$lambda
  beta0 <- result$beta[1, 1]
  alpha0 <- result$alpha
  K_gpr <- result$K
  u_weight <- result$u_hat
  base_est <- result$base_est
  
  noise_hat <- estimate_noise(Y, lambda, beta0, alpha0, K_gpr)
  
  sigma2_hat <- noise_hat$sigma2_hat
  tau_hat <- sigma2_hat / lambda
  
  test <- match.arg(test, c("asym", "boot"))
  func_name <- paste0("test_", test)
  
  pvalue <- do.call(func_name, list(n = n, Y = Y, X12 = X12, 
                                    beta0 = beta0, alpha0 = alpha0,
                                    K_gpr = K_gpr, sigma2_hat = sigma2_hat, 
                                    tau_hat = tau_hat, B = B))
  
  if (!is.null(data_test) & !is.null(fit_test)) {
    
    kern_size <- length(kern_list)
    Yhat_test <- 0
    for (d in seq(kern_size)) {
      kern <- kern_list[[d]]
      K1_m <- kern(fit_test$X1, X1)
      K2_m <- kern(fit_test$X2, X2)
      K_test <- K1_m + K2_m
      
      Yhat_base <- base_est$beta_list[[d]][1, 1] + 
        K_test %*% base_est$alpha_list[[d]]
      Yhat_test <- Yhat_test + u_weight[d] * Yhat_base
    }
    
    SSE_test <- sum((fit_test$Y - Yhat_test) ^ 2)
    test_RMSE <- sqrt(SSE_test / n)
    
  } else {
    test_RMSE <- -1
  }
  
  SSE_train <- sum((Y - beta0 - K_gpr %*% alpha0) ^ 2)
  train_RMSE <- sqrt(SSE_train / n)
  
  V0_inv <-
    compute_stat(n, Y, X12, beta0, K_gpr, sigma2_hat, tau_hat)$V0_inv
  
  K_eig <- svd(K_gpr)$d
  
  list(pvalue = pvalue, u_weight = u_weight, 
       lambda = lambda, train_RMSE = train_RMSE, 
       test_RMSE  = test_RMSE, trtst_ratio = train_RMSE / test_RMSE, 
       K_tr = tr(K_gpr), A_tr = tr(noise_hat$A), 
       V0_inv_tr = tr(V0_inv), K_eig = K_eig)
}



#' Conducting Score Tests for Interaction Using Asymptotic Test
#'
#' Conduct score tests comparing a fitted model and a more general alternative
#' model using asymptotic test.
#'
#' \bold{Asymptotic Test}
#'
#' This is based on the classical variance component test to construct a
#' testing procedure for the hypothesis about Gaussian process function.
#'
#' @param n (integer) A numeric number specifying the number of observations.
#' @param Y (vector of length n) Reponses of the dataframe.
#' @param X12 (dataframe, n*(p1\*p2)) The interaction items of first and second 
#' types of factors in the dataframe.
#' @param beta0 (numeric) Estimated bias of the model.
#' @param alpha0 (vector of length n) Estimated coefficients of the estimated 
#' ensemble kernel matrix.
#' @param K_gpr (matrix, n*n) Estimated ensemble kernel matrix.
#' @param sigma2_hat (numeric) The estimated noise of the fixed effects.
#' @param tau_hat (numeric) The estimated noise of the random effects.
#' @param B (integer) A numeric value indicating times of resampling 
#' when test = "boot".
#' @return \item{pvalue}{(numeric) p-value of the test.}
#' @author Wenying Deng
#' @seealso method: \code{\link{generate_kernel}}
#'
#' mode: \code{\link{tuning}}
#'
#' strategy: \code{\link{ensemble}}
#' @references Xihong Lin. Variance component testing in generalised linear
#' models with random effects. June 1997.
#'
#' Arnab Maity and Xihong Lin. Powerful tests for detecting a gene effect in
#' the presence of possible gene-gene interactions using garrote kernel
#' machines. December 2011.
#'
#' Petra Bu ̊zˇkova ́, Thomas Lumley, and Kenneth Rice. Permutation and
#' parametric bootstrap tests for gene-gene and gene-environment interactions.
#' January 2011.
#' 
#' @export test_asym
test_asym <- function(n, Y, X12, beta0, alpha0,
                      K_gpr, sigma2_hat, tau_hat, B) {
  
  score_chi <-
    compute_stat(n, Y, X12, beta0, K_gpr, sigma2_hat, tau_hat)$test_stat
  
  K0 <- K_gpr
  K12 <- X12 %*% t(X12)
  V0_inv <- ginv(tau_hat * K0 + sigma2_hat * diag(n))
  one <- rep(1, n)
  P0_mat <- V0_inv - V0_inv %*%
    one %*% ginv(t(one) %*% V0_inv %*% one) %*% t(one) %*% V0_inv
  
  drV0_tau <- K0
  drV0_sigma2 <- diag(n)
  drV0_del <- tau_hat * K12
  
  I0 <- compute_info(P0_mat,
                     mat_del = drV0_del, mat_sigma2 = drV0_sigma2,
                     mat_tau = drV0_tau)
  
  #Effective Info for delta
  tot_dim <- ncol(I0)
  I_deldel <-
    I0[1, 1] -
    I0[1, 2:tot_dim] %*% ginv(I0[2:tot_dim, 2:tot_dim]) %*% I0[2:tot_dim, 1]
  
  md <- tau_hat * tr(K12 %*% P0_mat) / 2
  
  m_chi <- I_deldel / (2 * md)
  d_chi <- md / m_chi
  
  pvalue <- 1 - pchisq(score_chi / m_chi, d_chi)
  
  pvalue
}



#' Conducting Score Tests for Interaction Using Bootstrap Test
#'
#' Conduct score tests comparing a fitted model and a more general alternative
#' model using bootstrap test.
#'
#' \bold{Bootstrap Test}
#'
#' When it comes to small sample size, we can use bootstrap test instead, which
#' can give valid tests with moderate sample sizes and requires similar
#' computational effort to a permutation test.
#'
#' @param n (integer) A numeric number specifying the number of observations.
#' @param Y (vector of length n) Reponses of the dataframe.
#' @param X12 (dataframe, n*(p1\*p2)) The interaction items of first and second 
#' types of factors in the dataframe.
#' @param beta0 (numeric) Estimated bias of the model.
#' @param alpha0 (vector of length n) Estimated coefficients of the estimated 
#' ensemble kernel matrix.
#' @param K_gpr (matrix, n*n) Estimated ensemble kernel matrix.
#' @param sigma2_hat (numeric) The estimated noise of the fixed effects.
#' @param tau_hat (numeric) The estimated noise of the random effects.
#' @param B (integer) A numeric value indicating times of resampling 
#' when test = "boot".
#' @return \item{pvalue}{(numeric) p-value of the test.}
#' @author Wenying Deng
#' @seealso method: \code{\link{generate_kernel}}
#'
#' mode: \code{\link{tuning}}
#'
#' strategy: \code{\link{ensemble}}
#' @references Xihong Lin. Variance component testing in generalised linear
#' models with random effects. June 1997.
#'
#' Arnab Maity and Xihong Lin. Powerful tests for detecting a gene effect in
#' the presence of possible gene-gene interactions using garrote kernel
#' machines. December 2011.
#'
#' Petra Bu ̊zˇkova ́, Thomas Lumley, and Kenneth Rice. Permutation and
#' parametric bootstrap tests for gene-gene and gene-environment interactions.
#' January 2011.
#' 
#' @export test_boot
test_boot <- function(n, Y, X12, beta0, alpha0,
                      K_gpr, sigma2_hat, tau_hat, B) {
  
  meanY <- K_gpr %*% alpha0 + beta0
  # lam_star <- sigma2_hat / tau_hat
  bs_test <- sapply(1:B, function(k) {
    
    Ystar <- meanY + rnorm(n, sd = sqrt(sigma2_hat))
    compute_stat(n, Ystar, X12, beta0, K_gpr, sigma2_hat, tau_hat)$test_stat
  })
  
  # assemble test statistic
  original_test <-
    compute_stat(n, Y, X12, beta0, K_gpr, sigma2_hat, tau_hat)$test_stat
  
  pvalue <- sum(as.numeric(original_test) <= bs_test) / B
  
  pvalue
}

