\name{stls}
\alias{stls}
\alias{print,stls-method}
\alias{summary,stls-method}
\alias{print,summary.stls-method}
\alias{coef,stls-method}
\alias{vcov,stls-method}
\alias{residuals,stls-method}
\alias{fitted,stls-method}

\title{
Estimation of truncated regression models using the Symmetrically Trimmed Least Squares (STLS) estimator
}
\description{
Function for estimation of linear regression models with truncated response variables (fixed truncation point), using the STLS estimator (Powell 1986)
}
\usage{
stls(formula, data, point = 0, direction = "left", beta = "ols", 
   covar = FALSE, na.action, ...)
\S4method{print}{stls}(x, digits = max(3, getOption("digits") - 2),
    width= getOption("width"), ...)
\S4method{summary}{stls}(object, level=0.95, ...)
\S4method{print}{summary.stls}(x, digits= max(3, getOption("digits") - 2), 
    width= getOption("width"), ...)
\S4method{coef}{stls}(object,...)
\S4method{vcov}{stls}(object,...)
\S4method{residuals}{stls}(object,...)
\S4method{fitted}{stls}(object,...)
}
\arguments{
\item{x, object}{
     an object of class \code{"stls"}
}
  \item{formula}{
     a symbolic description of the model to be estimated
}
  \item{data}{
     an optional data frame
}
  \item{point}{
     the value of truncation (the default is 0)
}
  \item{direction}{
     the direction of truncation, either \code{"left"} (the default) or \code{"right"}
}
 \item{beta}{the method of determining the starting values of the regression coefficients (See Details for more information):
\itemize{
     \item The default method is \code{"ols"}, meaning that the estimated regression coefficients from fitting a linear model with \code{\link{lm}} are used. 
\item Method \code{"ml"} means that the estimated regression coefficients from fitting a maximum likelihood model for truncated regression, assuming gaussian errors, are used. The maximum likelihood model is fitted using \code{\link[truncreg]{truncreg}}. 
\item The third option is to manually provide starting values as either a vector, column matrix or row matrix. 
}
}
  \item{covar}{
logical. Indicates whether the covariance matrix should be estimated. If \code{TRUE} the covariance matrix is estimated using bootstrap. The default number of replicates is 200 but this can be adjusted (see argument \code{\dots}). However, since the bootstrap procedure is time-consuming the default is \code{covar=FALSE}.
}
  \item{na.action}{
a function which indicates what should happen when the data contain \code{NA}s.
}
 \item{digits}{
    the number of digits to be printed     
}
 \item{width}{
     the width of the printing
}
 \item{level}{
     the desired level of confidence, for confidence intervals provided by \code{summary.stls}. A number between 0 and 1. The default value is \code{0.95}.
}
  \item{\dots}{
additional arguments. For \code{stls} the number of bootstrap replicates can be adjusted by setting \code{R=}the desired number of replicates. Also the \code{control} argument of \code{\link{optim}} can be set by \code{control=list()}.
}
}
\details{
Uses \code{\link{optim}} ("Nelder--Mead" method) to minimize the objective function described in Powell (1986) wrt the vector of regression coefficients in order to find the STLS estimates. The maximum number of iterations is set at 2000, but this can be adjusted by setting \code{control=list(maxit=\dots)} (see the \dots--argument). \cr\cr As the starting values of the regression coefficients have a great impact on the result of the minimization it is recommended to use one of the methods for generating these rather than supplying the values manually (unless one is confident that one has a good idea of what the starting values should be).
}
\value{
\code{stls} returns an object of class \code{"stls"}. \cr\cr The function \code{summary} prints a summary of the results, including two types of confidence intervals (t-distribution and bootstrap). The generic accessor functions
\code{coef}, \code{fitted}, \code{residuals} and \code{vcov} extract various useful features of the value returned by \code{stls}\cr\cr An object of class \code{"stls"}, a list with elements:
  \item{coefficients }{the named vector of coefficients}
  \item{startcoef }{the starting values of the regression coefficients used by \code{\link{optim}}}
  \item{value }{the value of the objective function corresponding to \code{coefficients}}
  \item{counts }{number of iterations used by \code{\link{optim}}. See the documentation for \code{\link{optim}} for further details}
  \item{convergence }{from \code{\link{optim}}. An integer code. 0 indicates successful completion. Possible error codes are \cr 1  indicating that the iteration limit maxit had been reached.\cr 10  indicating degeneracy of the Nelder--Mead simplex.}	
  \item{message }{from \code{\link{optim}}. A character string giving any additional information returned by the optimizer, or \code{NULL}.}
  \item{residuals }{the residuals of the model}
  \item{fitted.values }{the fitted values}
  \item{df.residual }{the residual degrees of freedom}
  \item{call }{the matched call}
  \item{covariance }{if \code{covar}\code{=}\code{TRUE}, the estimated covariance matrix}
  \item{bootrepl }{if \code{covar}\code{=}\code{TRUE}, the bootstrap replicates}
}
\references{
Powell, J. (1986) Symmetrically Trimmed Least Squares Estimation for Tobit Models, \emph{Econometrika}, \bold{54(6)}, pp 1435--1460
}
\author{
Anita Lindmark  and Maria Karlsson
}
\seealso{
\code{\link{stls.fit}}, the function that does the actual fitting \cr\cr
\code{\link{qme}}, for estimation of models with truncated response variables using the QME estimator \cr\cr
\code{\link{lt}}, for estimation of models with truncated response variables using the LT estimator \cr\cr
\code{\link[truncreg]{truncreg}} for estimating models with truncated response variables by maximum likelihood, assuming gaussian errors
}
\examples{
##Simulate a data.frame
n <- 10000
x1 <- runif(n,0,10)
x2 <- runif(n,0,10)
x3 <- runif(n,-5,5)
x4 <- runif(n,5,10)
x5 <- runif(n,-5,5)
y <- 1-2*x1+x2+2*x3+x4-x5+rnorm(n,0,2)
d <- data.frame(y=y,x1=x1,x2=x2,x3=x3,x4=x4,x5=x5)

##Use a truncated subsample
dtrunc <- subset(d, y>0)
  
##Use stls to estimate the model
stls(y~x1+x2+x3+x4+x5, dtrunc, point=0, direction="left", beta="ols", covar=FALSE)


##Example using data "PM10trunc"
data(PM10trunc)

stlspm10 <- 
stls(PM10~cars+temp+wind.speed+temp.diff+wind.dir+hour+day, data=PM10trunc, point=2)

summary(stlspm10)

}

\keyword{ regression }
