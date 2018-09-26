#' Conducting Gaussian Process Regression
#'
#' Conduct gaussian process regression based on the estimated ensemble kernel
#' matrix.
#'
#' After obtaining the ensemble kernel matrix, we can calculate the outpur of
#' gaussian process regression, the solution is given by
#' \deqn{\hat{\beta}=[1^T(K+\lambda I)^{-1}1]^{-1}1^T(K+\lambda I)^{-1}y}
#' \deqn{\hat{\alpha}=(K+\lambda I)^{-1}(y-\hat{\beta}1)} where
#' \eqn{\beta=intercept}.
#'
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
#' @param lambda (numeric) A numeric string specifying the range of 
#' noise to be chosen. The lower limit of lambda must be above 0.
#' @return \item{lam}{(numeric) The selected tuning parameter based on the 
#' estimated ensemble kernel matrix.}
#'
#' \item{intercept}{(numeric) Estimated bias of the model.}
#'
#' \item{alpha}{(vector of length n) Estimated coefficients of the estimated 
#' ensemble kernel matrix.}
#'
#' \item{K}{(matrix, n*n) Estimated ensemble kernel matrix.}
#'
#' \item{u_hat}{(vector of length K) A vector of weights of the kernels in 
#' the library.}
#' @author Wenying Deng
#' @seealso strategy: \code{\link{ensemble}}
#' @examples
#'
#'
#' estimation(Y, X1, X2, kern_list, mode = "loocv", strategy = "erm",
#' beta = 1, lambda = exp(seq(-5, 5)))
#'
#'
#' @export estimation
estimation <- function(
  Y, X1, X2, kern_list, mode = "loocv", strategy = "erm", 
  beta = 1, lambda_list = exp(seq(-10, 5))) {
  # base model estimate
  n <- length(Y)
  kern_size <- length(kern_list)
  base_est <- estimate_base(n, kern_size, Y, X1, X2, kern_list, 
                            mode, lambda_list)
  
  # ensemble estimate
  P_K_hat <- base_est$P_K_hat
  error_mat <- base_est$error_mat
  ens_res <- ensemble(n, kern_size, strategy, beta, error_mat, P_K_hat)
  
  # assemble ensemble kernel matrix
  K_ens <- ensemble_kernel_matrix(ens_res$A_est)
  
  # final estimate
  lambda_ens <- tuning(Y, K_ens, mode, lambda_list)
  ens_est <- estimate_ridge(X = matrix(1, nrow = n, ncol = 1),
                            K = K_ens, Y = Y, lambda = lambda_ens)
  list(lambda = lambda_ens, 
       beta = ens_est$beta, 
       alpha = ens_est$alpha, K = K_ens, 
       u_hat = ens_res$u_hat, base_est = base_est)
}

#' Estimating Projection Matrices for all base models
#' 
#' Calculate the estiamted projection matrices for every kernels in the kernel
#' library.
#' 
#' For a given mode, this function return a list of projection matrices for
#' every kernels in the kernel library and a n*kern_size matrix indicating
#' errors.
#' 
#' @param n (integer) A numeric number specifying the number of observations.
#' @param kern_size (integer, =K) A numeric number specifying the number of 
#' kernels in the kernel library.
#' @param Y (vector of length n) Reponses of the dataframe.
#' @param X1 (dataframe, n*p1) The first type of factor in the dataframe (could 
#' contains several subfactors).
#' @param X2 (dataframe, n*p2) The second type of factor in the dataframe (could 
#' contains several subfactors).
#' @param kern_list (list of length K) A list of kernel functions given by user.
#' @param mode (character) A character string indicating which tuning parameter 
#' criteria is to be used.
#' @param lambda (numeric) A numeric string specifying the range of noise 
#' to be chosen. The lower limit of lambda must be above 0.
#' @return \item{A_hat}{(list of length K) A list of projection matrices for 
#' every kernels in the kernel library.}
#' 
#' \item{error_mat}{(matrix, n*K) A n\*kern_size matrix indicating errors.}
#' @author Wenying Deng
#' @references Jeremiah Zhe Liu and Brent Coull. Robust Hypothesis Test for
#' Nonlinear Effect with Gaus- sian Processes. October 2017.
#' @examples
#' 
#' 
#' estimate_base(n = 100, kern_size = 3, Y, X1, X2, kern_list,
#' mode = "loocv", lambda = exp(seq(-5, 5)))
#' 
#' 
#' @export estimate_base
estimate_base <- function(n, kern_size, Y, X1, X2, kern_list, mode, lambda){
  A_hat <- list()
  P_K_hat <- list()
  beta_list <- list()
  alpha_list <- list()
  lambda_list <- list()
  error_mat <- matrix(0, nrow = n, ncol = kern_size)
  
  for (d in seq(kern_size)) {
    kern <- kern_list[[d]]
    K1_m <- kern(X1, X1)
    K2_m <- kern(X2, X2)
    K <- K1_m + K2_m
    K_scale <- tr(K)
    K <- K / K_scale
    
    if (length(lambda) != 0) {
      lambda0 <- tuning(Y, K, mode, lambda)
      estimate <- estimate_ridge(X = matrix(1, nrow = n, ncol = 1),
                                 K = K, Y = Y, lambda = lambda0)
      A <- estimate$proj_matrix$total
      
      # produce loocv error matrix
      error_mat[, d] <- (diag(n) - A) %*% Y / (1 - diag(A))
      
      A_hat[[d]] <- A
      P_K_hat[[d]] <- estimate$proj_matrix$P_K0
      beta_list[[d]] <- estimate$beta
      alpha_list[[d]] <- estimate$alpha
      lambda_list[[d]] <- lambda0
    }
  }
  
  list(A_hat = A_hat, P_K_hat = P_K_hat,
       beta_list = beta_list, alpha_list = alpha_list,
       lambda_list = lambda_list, error_mat = error_mat)
}


#' Estimating Projection Matrices and Parameter Estimates for an Single Model
estimate_ridge <- function(X, K, Y, lambda){
  # standardize kernel matrix
  n <- nrow(K)
  
  V_inv <- ginv(K + lambda * diag(n))
  B_mat <- ginv(t(X) %*% V_inv %*% X) %*% t(X) %*% V_inv
  
  # project matrices
  P_X <- X %*% B_mat  # projection to fixed-effect space
  P_K0 <- K %*% V_inv # projection to kernel space
  P_K <- P_K0 %*% (diag(n) - P_X) # residual projection to kernel space
  
  # parameter estimates
  beta <- B_mat %*% Y
  alpha <- V_inv %*% (diag(n) - P_X) %*% Y
  
  # return
  proj_matrix_list <- list(total = P_X + P_K, 
                           P_X = P_X, P_K = P_K, P_K0 = P_K0)
  list(beta = beta, alpha = alpha, 
       proj_matrix = proj_matrix_list)
}

#' Estimating Ensemble Kernel Matrix
ensemble_kernel_matrix <- function(A_est, eig_thres = 1e-9){
  As <- svd(A_est)
  U <- As$u
  d <- As$d
  
  # produce spectral components for ensemble kernel matrix
  ensemble_dim <- sum(d > eig_thres)
  U_ens <- U[, 1:ensemble_dim]
  d_ens <- d[1:ensemble_dim]/(1 - d[1:ensemble_dim])
  
  # assemble ensemble matrix and return
  K_hat <- U_ens %*% diag(d_ens) %*% t(U_ens)
  K_hat / tr(K_hat)
}
