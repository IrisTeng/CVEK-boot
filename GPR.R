kernelReg <- 
  function(formula, data = NULL, X_predict=NULL,label_names = label_names,l=1,sigma2_f=1,sigma2_n=1,int=FALSE){
    # prepare data
    y <- data[,as.character(attr(terms(formula),"variables"))[2]]
    
    # extract the information from the given formula to construct a true formula
    generic_formula0 <- 
      generic_formula(formula=formula,label_names = label_names)$generic_formula
    
    # extract design matrix
    X <- model.matrix(generic_formula0 ,data=data)[,-1]
    
    if(int==TRUE){
      X1 <- X[,1:length(label_names[[1]])]
      X2 <- X[,(length(label_names[[1]])+1):dim(X)[2]]
      X_int <- NULL
      for (i in 1:length(label_names[[1]])){
        X_int <- cbind(X_int, X1[,i]*X2)
      }
      X <- X_int
    }
    
    # estimation
    SE <- function(xp,xq,l,sigma2_f) sigma2_f*exp(-sum((xp-xq)^2)/(2*l^2))
    
    cov <- function(X1,X2) apply(X1,1,function(xp){
      apply(X2,1,function(xq){
        SE(xp,xq,l=l,sigma2_f=sigma2_f)
      })
    })

    # use equation (2.35)
    K <- cov(X,X)
    A <- solve(K+diag(dim(X)[1])*sigma2_n)
    f_bar <- K%*%A%*%y
    
    # prediction, equations (2.25) and (2.26)
    if(!is.null(X_predict)){
      K_star <- cov(X_predict,X)
      f_predict <- t(K_star)%*%A%*%y
      cov_predict <- cov(X_predict,X_predict)-t(K_star)%*%A%*%K_star
    }
    else{
      K_star <- NULL
      f_predict <- NULL
      cov_predict <- NULL
    }
    return(list(f_bar=f_bar,f_predict=f_predict,cov_predict=cov_predict))
  }


# compare with gausspr
x_predict <- matrix(rnorm(20),5,4)
label_names <- list(X1=c("x1","x2"),X2=c("x3","x4"))
rawData <- originalData(2,5,1,3,20,label_names=label_names,beta_int=0.1,scale=10,eps=1)
K0 <- kernelReg(Y~X1+X2,data=rawData,X_predict=x_predict,
                label_names=label_names,l=2,sigma2_f=1,sigma2_n=0.6)
K0$cov_predict

gp <- gausspr(rawData[,-1],rawData[,1], kernel="rbfdot", 
              kpar=list(sigma=1/8), fit=T,scaled=FALSE, var=0.6)
cbind(K0$f_bar,fitted(gp))
cbind(K0$f_predict,predict(gp,x_predict))
