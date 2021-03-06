% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernelGenerate.R
\name{kernelGenerate}
\alias{kernelGenerate}
\title{Generating A Single Kernel}
\usage{
kernelGenerate(method = "rbf", Sigma = 0, l = 1, p = 2)
}
\arguments{
\item{method}{A character string indicating which kernel is to be computed.}

\item{Sigma}{The covariance matrix for neural network kernel.}

\item{l}{A numeric number indicating the hyperparameter (flexibility) of a
specific kernel.}

\item{p}{For polynomial, p is the power; for matern, v = p + 1 / 2; for
rational, alpha = p.}
}
\value{
\item{Kern}{A function indicating the generated kernel.}
}
\description{
Generate kernels for the kernel library.
}
\details{
There are seven kinds of kernel available here. For convenience, we define
\eqn{r=\mid x-x'\mid}.

\bold{Gaussian RBF Kernels} \deqn{k_{SE}(r)=exp\Big(-\frac{r^2}{2l^2}\Big)}

\bold{Matern Kernels}
\deqn{k_{Matern}(r)=\frac{2^{1-\nu}}{\Gamma(\nu)}\Big(\frac{\sqrt{2\nu
r}}{l}\Big)^\nu K_\nu \Big(\frac{\sqrt{2\nu r}}{l}\Big)}

\bold{Rational Quadratic Kernels} \deqn{k_{RQ}(r)=\Big(1+\frac{r^2}{2\alpha
l^2}\Big)^{-\alpha}}

\bold{Polynomial Kernels} \deqn{k(x, x')=(x \cdot x')^p} We have intercept
kernel when \eqn{p=0}, and linear kernel when \eqn{p=1}.

\bold{Neural Network Kernels} \deqn{k_{NN}(x,
x')=\frac{2}{\pi}sin^{-1}\Big(\frac{2\tilde{x}^T\Sigma
\tilde{x}'}{\sqrt{(1+2\tilde{x}^T\Sigma \tilde{x})(1+2\tilde{x}'^T\Sigma
\tilde{x}')}}\Big)}
}
\examples{


##kernelGenerate(method = "rbf", Sigma = 0, l = 1, p = 2)

##Kernlist <- NULL
##Kernlist <- c(Kernlist, kernelGenerate('rbf', l = .6))
##Kernlist <- c(Kernlist, kernelGenerate('rbf', l = 1))
##Kernlist <- c(Kernlist, kernelGenerate('rbf', l = 2))


}
\references{
The MIT Press. Gaussian Processes for Machine Learning, 2006.
}
\author{
Wenying Deng
}
