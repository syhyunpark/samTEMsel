% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/samTEMsel.main.R
\name{samTEMsel}
\alias{samTEMsel}
\title{Sparse Additive Models for Treatment Effect-Modifier Selection (main function)}
\usage{
samTEMsel(y, A, X, mu.hat = NULL, d = 3, lambda = NULL,
  nlambda = 50, lambda.min.ratio = 0.01, thol = 1e-05,
  max.ite = 1e+05, regfunc = "L1", terms.fit = FALSE, basis = NULL,
  basisc = NULL)
}
\arguments{
\item{y}{a n-by-1 vector of responses}

\item{A}{a n-by-1 vector of treatment variable; each element represents one of the L(>1) treatment conditions; e.g., c(1,2,1,1,3...); can be a factor-valued}

\item{X}{a n-by-p matrix of pretreatment features}

\item{mu.hat}{a n-by-1 vector of the fitted X main effect term of the model provided by the user; defult is \code{NULL}, in which case \code{mu.hat} is taken to be a vector of zeros; the optimal choice for this vector is E(y|X)}

\item{d}{number of basis spline functions to be used for each component function; the default value is 3; d=1 corresponds to the linear model}

\item{lambda}{a user-supplied lambda sequence; typical usage is to have the program compute its own lambda sequence based on \code{nlambda} and \code{lambda.min.ratio}.}

\item{nlambda}{total number of \code{lambda} values; the default value is 50.}

\item{lambda.min.ratio}{the smallest value for lambda, as a fraction of lambda.max, the (data derived) entry value (i.e. the smallest value for which all coefficients are zero); the default is 0.01.}

\item{thol}{stopping precision; the default value is 1e-5.}

\item{max.ite}{number of maximum iterations; the default value is 1e5.}

\item{regfunc}{type of the regularizer; the default is "L1"; can also be "MCP" or "SCAD".}

\item{terms.fit}{if \code{TRUE}, return the component-wise fitted vectors (\code{y.hat.list}) and the partial residuals (\code{resid.list}), evaluated over a grid of \code{lambda}.}

\item{basis}{default is \code{NULL}; one can provide a n-by-d*p*L basis matrix associted with the spline representation of the model; this is to efficiently implement cross-validation in \code{cv.samTEMsel}.}

\item{basisc}{default is \code{NULL}; one can provide a n-by-d*p*(L-1) basis matrix associted with the spline representation that incorporates the "orthogonality" constraint; this is only to efficiently implement cross-validation in \code{cv.samTEMsel}.}
}
\value{
a list of information of the fitted constrained additive models including
 \item{w}{the solution path matrix; the estimated d*p*L-by-nlambda coefficient matrix associted with B-splines, with each column corresponding to a regularization parameter in \code{lambda}}
 \item{wc}{the solution path matrix; the estimated d*p*(L-1)-by-nlambda coefficient matrix associted with B-splines that incorports the "orthogonality" constraint, with each column corresponding to a regularization parameter in \code{lambda}}
 \item{lambda}{a sequence of regularization parameter}
 \item{d}{the number of basis spline functions used for each component function}
 \item{func_norm}{the functional norm matrix (p-by-nlambda) with each column corresponding to a regularization parameter in \code{lambda}; since we have p variables, the length of each column is p.}
 \item{df}{the degree of freedom of the solution path (the number of non-zero component functions).}
 \item{knots}{a (d-1)-by-p matrix; each column contains the knots applied to the corresponding variable.}
 \item{Boundary.knots}{a 2-by-p matrix; each column contains the boundary points applied to the corresponding variable.}
 \item{sse}{sums of square errors of the fitted models over the solution path.}
 \item{intercept}{the list of L treatment-specific intercepts; each element of the list is a 1-by-nlambda matrix, with each column corresponding to a regularization parameter in \code{lambda}.}
 \item{X.min}{a p-by-1 vector, with each entry corresponding to the minimum of each input variable; used for rescaling in testing.}
 \item{X.ran}{a p-by-1 vector, with each entry corresponding to the range of each input variable; used for rescaling in testing.}
 \item{residuals}{the n-by-nlambda matrix of residuals of the fitted models, with each column corresponding to each regularization parameter in \code{lambda}.}
 \item{y.hat}{the n-by-nlambda matrix of fitted values for y, with each column corresponding to each regularization parameter in \code{lambda}.}
 \item{resid.list}{the list of p component-wise partial residuals, in which each (the jth) element is a n-by-nlambda matrix of the jth partial residuals, with each column corresponding to each regularization parameter in \code{lambda}.}
 \item{y.hat.list}{the list of p component-wise fitted values, in which each (the jth) element is a n-by-nlambda matrix of the jth fitted values, with each column corresponding to each regularization parameter in \code{lambda}.}
 \item{basis}{the n-by-d*p*L basis matrix associted with the spline representation of the (unconstrained) additive model.}
 \item{basisc}{the n-by-d*p*(L-1) basis matrix associted with the spline representation of the model that incorporates the "orthogonality" constraint.}
}
\description{
The function \code{samTEMsel} implements estimation of a constrained sparse additve model.
}
\details{
A constrained additive regression model represents the joint effects of treatment and p pretreatment covariates on an outcome via treatment-specific additive component functions defined over the p covariates, subject to the constraint that the expected value of the outcome given the covariates equals zero, while leaving the main effects of the covariates unspecified.
Under this flexible representation, the treatment-by-covariates interaction effects are determined by distinct shapes (across treatment levels) of the unspecified component functions.
Optimized under a penalized least square criterion with a L1 (or SCAD/MCP) penalty, the constrained additive model can effectively identify/select treatment effect-modifiers (from the p pretreatment covariates) that exhibit possibly nonlinear interactions with the treatment variable; this is achieved by producing a sparse set of estimated component functions.
The estimated nonzero component functions (available from the returned \code{samTEMsel} object) can be used to make individualized treatment recommendations (ITRs) for future subjects; see also \code{make_ITR} for such ITRs.

The regularization path is computed at a grid of values for the regularization parameter \code{lambda}.
}
\examples{
set.seed(112)
n.train = 400
n.test  = 50
n = n.train + n.test
p = 20
A = rbinom(n, 1, 0.5) + 1  # treatment variable taking a value in {1,2} with equal prob.
X = matrix(runif(n*p, -pi/2,pi/2), n, p)  # pretreatment covariates
noise = rnorm(n, 0, 0.5)
# X main effect on y; a highly nonlinear (cosine) function; depends on the first 10 covariates
main.effect = rep(0, n); for(j in 1:10){
  main.effect = main.effect + cos(X[,j])
}
# A-by-X ineraction effect on y; depends only on X1 and X2.
interaction.effect = (A-1.5)*X[,1]  +  2 * (A-1.5)*(cos(X[,2]) - 0.5)
# generate outcome y
y = main.effect  + interaction.effect + noise


# train/test set splitting
train.index = 1:n.train
y.train = y[train.index]
X.train = X[train.index,]
A.train = A[train.index]
y.test = y[-train.index]
X.test = X[-train.index,]
A.test = A[-train.index]

# fit samTEMsel() based on the training set
samTEMsel.obj = samTEMsel(y.train, A.train, X.train, nlambda = 50)

# a n.test-by-nlambda matrix of predicted values:
predict_samTEMsel(samTEMsel.obj, newX = X.test, newA = A.test)

# pick a particular lambda.index, say, 10, as the regularization parameter.
lambda.index = 10
# for an optimal selection of lambda.index, see cv.samTEMsel().

# a n.test-by-1 vector of predicted values (given lambda.index = 10)
predict_samTEMsel(samTEMsel.obj, newX = X.test, newA = A.test, lambda.index=lambda.index)

# the estimated L2 norm of the component functions (given lambda.index=10)
samTEMsel.obj$func_norm[,lambda.index]

# p component-wise fitted values (given lambda.index=10); the last column is the intercept
predict_samTEMsel(samTEMsel.obj, X.test, A.test, type="terms", lambda.index =lambda.index)

# can plot the estimated component functions (say, the first two functions, j=1,2)
plot_samTEMsel(samTEMsel.obj, which.index = c(1,2), lambda.index = lambda.index)

# can make inidividualized treatment recommendations (ITRs)
trt.rule = make_ITR(samTEMsel.obj, newX = X.test, lambda.index = lambda.index)$trt.rule
head(trt.rule)

# an (IPWE) estimate of the "value" of this particualr treatment rule, trt.rule:
mean(y.test[A.test==trt.rule])

# compare the above value to the following estimated "values" of "naive" treatment rules:
mean(y.test[A.test==1])   # just assign everyone to A=1
mean(y.test[A.test==2])   # just assign everyone to A=2

}
\seealso{
\code{cv.samTEMsel}, \code{predict_samTEMsel}, \code{plot_samTEMsel}, \code{make_ITR}
}
\author{
Park, Petkova, Tarpey, Ogden
}
