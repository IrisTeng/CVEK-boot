# Hello World!
# from vectors to single variables
generic_formula <- 
  function(formula, label_names = label_names){
    formula_factors <- attr(terms(formula),"factors")
    generic_formula <- Y~1
    length_main_effect <- 0
    for(i in 1:dim(formula_factors)[2]){
      terms_names <- 
        rownames(formula_factors)[which(formula_factors[,i]==1)]
      if(length(terms_names)==1){
        generic_formula <- 
          update.formula(generic_formula,
                         as.formula(paste(as.character(
                           attr(terms(formula),"variables"))[2],
                           paste(label_names[[terms_names]], collapse=" + "), sep=" ~ .+")))
        length_main_effect <- length_main_effect+length(label_names[[terms_names]])
      }
      else{
        interaction_formula <- 
          paste("(",paste(label_names[[terms_names[1]]], collapse=" + "),")",sep="")
        for(j in 2:length(terms_names)){
          interaction_formula <- 
            paste(interaction_formula,"*(",
                  paste(label_names[[terms_names[j]]], collapse=" + "),")",sep=" ")
        }
        generic_formula <- 
          update.formula(generic_formula,
                         as.formula(paste(as.character(
                           attr(terms(formula),"variables"))[2], interaction_formula, sep=" ~ .+")
                         )
          )
      }
    }
    return(list(generic_formula=generic_formula,
                length_main_effect=length_main_effect))
  }

#generic_formula(Y~X1+X2+X3+X2*X3,list(X1=c("x1", "x2","x3"),X2=c("x4","x5"),X3=c("x6","x7")))

library(mvtnorm)
originalData2 <- 
  function(size, label_names = label_names, 
           method=NULL, int_effect=0, eps=0.1){
    
    # X~N(0,I)
    l <- 1
    X1 <- rmvnorm(n=size,
                  mean=rep(0,length(label_names[[1]])),
                  sigma=diag(length(label_names[[1]])))
    X2 <- rmvnorm(n=size,
                  mean=rep(0,length(label_names[[2]])),
                  sigma=diag(length(label_names[[2]])))
    
    # w~N(0,I)
    p <- length(label_names[[1]])+length(label_names[[2]])
    q <- length(label_names[[1]])*length(label_names[[2]])
    w0 <- rmvnorm(n=1, mean=rep(0,p),sigma=diag(p))
    
    if(int_effect!=0){
      X_int <- NULL
      for (i in 1:length(label_names[[1]])){
        X_int <- cbind(X_int, X1[,i]*X2)
      }
      X <- cbind(X1, X2, X_int)
      w <- cbind(w0, int_effect*rmvnorm(n=1, mean=rep(0,q),sigma=diag(q)))
    }
    else{
      X <- cbind(X1, X2)
      w <- w0
    }
    
    if(method=="linear")
      SE <- function(xp,xq,l) t(xp)%*%xq
    else if(method=="gaussian")
      SE <- function(xp,xq,l) exp(-sum((xp-xq)^2)/(2*l^2))
    else
      stop("method must be linear or gaussian")
    cov <- function(X1,X2) apply(X1,1,function(xp){
      apply(X2,1,function(xq){
        SE(xp,xq,l=l)
      })
    })

    Y <- rmvnorm(n=1, mean=rep(0,size),sigma=cov(X,X)+eps*diag(size))+rnorm(1)
    data <- as.data.frame(cbind(t(Y),X1,X2))
    colnames(data) <- c("Y", label_names[[1]],label_names[[2]])
    return(data)
  }

#label_names <- list(X1=c("x1","x2"),X2=c("x3","x4"))
# rawData <- originalData2(20,label_names,method = "linear", int_effect = 3, eps = 3)


# the lower limit of lambda must be above 0
library(MASS)
kernelReg <- 
  function(formula, data = NULL, label_names = label_names, 
           method=NULL, l=1, lambda=0){
    # prepare data
    y <- data[,as.character(attr(terms(formula),"variables"))[2]]
    
    # extract the information from the given formula to construct a true formula
    generic_formula0 <- 
      generic_formula(formula=formula,label_names = label_names)$generic_formula
    
    # extract design matrix
    X <- model.matrix(generic_formula0, data=data)[,-1]
    n <- nrow(X)
    
    # estimation
    if(method=="linear")
      SE <- function(xp,xq,l) t(xp)%*%xq
    else if(method=="gaussian")
      SE <- function(xp,xq,l) exp(-sum((xp-xq)^2)/(2*l^2))
    else
      stop("method must be linear or gaussian")
    cov <- function(X1,X2) apply(X1,1,function(xp){
      apply(X2,1,function(xq){
        SE(xp,xq,l=l)
      })
    })
    K <- cov(X,X)
    Loocv <- sapply(lambda, function(k){
      A <- K%*%ginv(K+k*diag(n))
      sum(((diag(n)-A)%*%y/diag(diag(n)-A))^2)
    })
    lambda0 <- lambda[which(Loocv==min(Loocv))]
    
    # estimate intercept
    K1 <- cbind(1,K)
    K2 <- cbind(0,rbind(0,K))
    # beta0 <- (solve(lambda0*K2+t(K1)%*%K1)%*%t(K1)%*%y)[1]
    
    # K is singular under linear case, so apply generalized inverse here
    theta <- ginv(lambda0*K2+t(K1)%*%K1)%*%t(K1)%*%y
    beta0 <- theta[1]
    alpha <- theta[-1]
    # use equation (2.35)
    # B <- K%*%solve(K+lambda0*diag(n))
    return(list(sigma2_n=lambda0, intercept=beta0, alpha=alpha))
  }

# kernelReg(Y~X1+X2, data = rawData, label_names, method = "linear",lambda  = seq(.1, 10, .1))


multiScore2 <- 
  function(formula, data = NULL, label_names = label_names, 
           lambda=0, method=NULL, l=1){
    
    if(method=="linear")
      SE <- function(xp,xq,l) t(xp)%*%xq
    else if(method=="gaussian")
      SE <- function(xp,xq,l) exp(-sum((xp-xq)^2)/(2*l^2))
    else
      stop("method must be linear or gaussian")
    cov <- function(X1,X2) apply(X1,1,function(xp){
      apply(X2,1,function(xq){
        SE(xp,xq,l=l)
      })
    })
    # prepare data and estimate under null
    y <- data[,as.character(attr(terms(formula),"variables"))[2]]
    re=generic_formula(formula=formula,label_names = label_names)
    generic_formula0 <- re[[1]]
    len <- re[[2]]
    X <- model.matrix(generic_formula0, data=data)[,-1]
    n <- nrow(X)
    X1 <- X[,c(1:length(label_names[[1]]))]
    X2 <- X[,c((length(label_names[[1]])+1):len)]
    X12 <- X[,c((len+1):dim(X)[2])]
    
    result <- kernelReg(Y~X1+X2, data, label_names, 
                        method, l, lambda)
    sigma2_n <- result[[1]]
    beta0 <- result[[2]]
    
    # compute score statistic
    K1 <- cov(X1,X1)
    K2 <- cov(X2,X2)
    K0 <- K1+K2
    K12 <- X12%*%t(X12)
    V0 <- K0+sigma2_n*diag(n)
    score <- t(y-beta0)%*%ginv(V0)%*%K12%*%ginv(V0)%*%(y-beta0)
    ret <- list(test.stat=score)
    return(ret)
  }
#multiScore2(Y~X1+X2+X1*X2, data = rawData, label_names, method = "linear",lambda  = seq(.1, 10, .1))

kernelBoot <- 
  function(formula, data = NULL, label_names = label_names, 
           method = NULL, B=100, l=1, lambda=seq(.1, 10, .01)){
    
    if(method=="linear")
      SE <- function(xp,xq,l) t(xp)%*%xq
    else if(method=="gaussian")
      SE <- function(xp,xq,l) exp(-sum((xp-xq)^2)/(2*l^2))
    else
      stop("method must be linear or gaussian")
    cov <- function(X1,X2) apply(X1,1,function(xp){
      apply(X2,1,function(xq){
        SE(xp,xq,l=l)
      })
    })
    result <- kernelReg(formula, data, label_names, 
                        method, l, lambda)
    sigma2_n <- result[[1]]
    beta0 <- result[[2]]
    alpha <- result[[3]]
    
    # prepare data
    y <- data[,as.character(attr(terms(formula),"variables"))[2]]
    
    # extract the information from the given formula to construct a true formula
    generic_formula0 <- 
      generic_formula(formula=formula,label_names = label_names)$generic_formula
    
    # extract design matrix
    X <- model.matrix(generic_formula0, data=data)[,-1]
    n <- nrow(X)
    inds <- 1:n
    
    # conduct bootstrap
    bs_test <- sapply(1:B, function(k){
      bs_ind <- sample(inds, size = n, replace = T)
      X_star <- data[bs_ind,-1]
      mean_Y <- cov(X_star,X)%*%alpha+beta0
      # cov_Y <- cov(X_star,X_star)-cov(X_star,X)%*%
      # ginv(cov(X,X)+sigma2_n*diag(n))%*%cov(X,X_star)
      Y <- rnorm(n, mean = mean_Y, sd = sqrt(sigma2_n))
      dat <- cbind(Y, X_star)
      multiScore2(Y~X1+X2+X1*X2,data=dat,
                  label_names, lambda, method)$test.stat
    })
    
    # assemble test statistic
    original_test <- multiScore2(Y~X1+X2+X1*X2,data=data,
                                 label_names, lambda, method)$test.stat

    bs_pvalue <- sum(as.numeric(original_test)<=bs_test)/B
    return(bs_pvalue)
  }


verify <- 
  function(i){
    rawData <- originalData2(size=n,label_names, method=method,int_effect = int_effect)
    linear_pvalue <- kernelBoot(formula, data = rawData, 
                                label_names = label_names, method = "linear", B=B, l=l, lambda=seq(.1, 10, .01))
    gaussian_pvalue <- kernelBoot(formula, data = rawData, 
                              label_names = label_names, method = "gaussian", B=B, l=l, lambda=seq(.1, 10, .01))
    # })
    # prob <- sum(pvalue<0.05)/M
    # return(prob)
    return(c(linear_pvalue,gaussian_pvalue))
  }

# produce plot similar to the CVEK paper. 
# (x-axis: interaction effect size, y-axis: ratio of p-value smaller than 0.05). 
# vary the interaction effect size between seq(0, 0.1, 0.02), 
# vary the effect size multiplication factor within  (1, 3, 5, 10), 
# vary the noise level within (0.05, 0.1, 0.25, 0.5, 1). 
# Keep bootstrap sample size at 200. Keep b11 = 2 (lower bound of effect size).
# number of repetition M fix at 1000.
# Repeat such experiment for linear and ridge regression.
library(snowfall)
library(mvtnorm)
library(MASS)
sfInit(parallel=T,cpus=20)
sfLibrary(mvtnorm)
sfLibrary(MASS)

formula <- Y~X1+X2
label_names <- list(X1=c("x1","x2"),X2=c("x3","x4"))
result <- NULL
for (M in c(1000)){
  for (method in c("linear","gaussian")){
    for(int_effect in seq(0,10,1)){
      for(l in c(0.3,0.6,1,3,5)){
        for(n in seq(200,1000,200)){
          for (B in c(200)){
            sfExport("formula","label_names","int_effect",
                     "n","B","M","l","method","kernelBoot","originalData2",
                     "generic_formula","kernelReg","multiScore2") 
            system.time(res <- sfSapply(1:M,verify))
            write.table(t(res),file="simulation_power_pvalue.txt",
                        row.names=F,col.names=F,append=T)            
            res2 <- apply(res,1,function(x){sum(x<0.05)/M})
            result <- rbind(result,c(int_effect,l,method,n,res2))
            cat(c(int_effect,l,method,n,res2),
                file="simulation_power.txt",append=T,"\n")
            cat("Finished:interaction effect size=",int_effect,"n=",n,"l=",l,"\n")
          }
        }
      }
    }
  }
}
write.csv(result,file="simulation_power.csv",row.names = F,quote = F)

# plot
# file <- read.csv(file="/Users/iristeng/Desktop/simulation_power.csv")
# #colnames(file) <- c('beta_int','b11','b12','b21','b22','n','B','prob')
# file$effect_ind <- ifelse(file$b11==2,1,
#                           ifelse(file$b11==6,2,ifelse(file$b11==10,3,4)))
# library(ggplot2)
# file$effect_ind <- factor(file$effect_ind, levels = 1:4,
#                           labels = c('effect size1','effect size2',
#                                      'effect size3','effect size4'))
# file$e <- factor(file$e, labels = c('e=0.05','e=0.10','e=0.25','e=0.50','e=1.00'))
# file$n <- as.factor(file$n)
# p <- ggplot(data = file, aes(x=beta_int, y=linear, shape=n, color=n))+
#   geom_point(size=0.8)+geom_line()+
#   geom_hline(yintercept=0.05)+facet_grid(effect_ind~e)+
#   labs(x='beta interaction', y='probability')+
#   theme_set(theme_bw())+theme(panel.grid=element_blank())
# p
# 
# q <- ggplot(data = file, aes(x=beta_int, y=ridge, shape=n, color=n))+
#   geom_point(size=0.8)+geom_line()+
#   geom_hline(yintercept=0.05)+facet_grid(effect_ind~e)+
#   labs(x='beta interaction', y='probability')+
#   theme_set(theme_bw())+theme(panel.grid=element_blank())
# q
