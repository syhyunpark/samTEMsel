#--------------------
# Package: samTEMsel
#--------------------


#' Sparse Additive Models for Treatment Effect-Modifier Selection (cross-validation function)
#'
#' Does k-fold cross-validation for \code{\link{samTEMsel}}, produces a plot, and returns the sequence of the fitted constrained additive models implied by the sequence of regularization parameters \code{lambda} and the index, \code{lambda.opt.index}, corresponding to the estimated optimal regularization parameter.
#'
#' @param y   a n-by-1 vector of responses
#' @param A   a n-by-1 vector of treatment variable; each element represents one of the L(>1) treatment conditions; e.g., c(1,2,1,1,3...); can be a factor-valued
#' @param X   a n-by-p matrix of pretreatment features
#' @param mu.hat   a n-by-1 vector of the fitted X main effect term of the model provided by the user; defult is \code{NULL}, in which case \code{mu.hat} is taken to be a vector of zeros; the optimal choice for this vector is E(y|X)
#' @param d     number of basis spline functions to be used for each component function; the default value is 3; d=1 corresponds to the linear model
#' @param n.folds  number of folds for cross-validation; the default is 10.
#' @param thol  stopping precision for the coordinate-descent algorithm; the default value is 1e-5.
#' @param max.ite  number of maximum iterations; the default value is 1e5.
#' @param regfunc  type of the regularizer; the default is "L1"; can also be "MCP" or "SCAD".
#' @param lambda   a user-supplied regularization parameter sequence; typical usage is to have the program compute its own lambda sequence based on nlambda and lambda.min.ratio.
#' @param nlambda  total number of lambda values; the default value is 50.
#' @param lambda.min.ratio  the smallest value for lambda, as a fraction of lambda.max, the (data derived) entry value (i.e. the smallest value for which all coefficients are zero); the default is 0.01.
#' @param row.ordering  a particular ordering (n-by-1 vector) of the rows (provided by the user) to be used for cross-validation; the default is \code{NULL} in which case the row order is randomly shuffled.
#' @param cv1sd  if \code{TRUE}, an optimal regularization parameter is chosen based on: the mean cross-validated error + 1 SD of the mean cross-validated error, which typically results in an increase in regularization; the defualt is \code{FALSE}.
#' @param terms.fit   if \code{TRUE}, in \code{samTEMsel}, compute and store the component-wise fitted vectors (y.hat.list) and partial residuals (resid.list), which are evaluated over a grid of \code{lambda}.
#' @param plots  if \code{TRUE}, produce a cross-validation plot of the estimated mean squared error versus the regulariation parameter index.
#'
#' @return a list of information of the fitted constrained sparse additive model including
#'  \item{samTEMsel.obj}{an object of class \code{samTEMsel}, which contains the sequence of the fitted constrained additive models implied by the sequence of the regularization parameters \code{lambda}; see \code{\link{samTEMsel}} for detail.}
#'  \item{lambda.opt.index}{an index number, indicating the index of the estimated optimal regularization parameter in \code{lambda}.}
#'  \item{nonzero.index}{a set of numbers, indicating the indices of estimated nonzero component functions, evalated at the regularization parameter index \code{lambda.opt.index}.}
#'  \item{func_norm.opt}{a p-by-1 vector, indicating the norms of the estimated component functions evaluatd at the regularization parameter index \code{lambda.opt.index}, with each element corresponding to the norm of each estimated component function.}
#'  \item{cv.storage}{a n.folds-by-nlambda matrix of the estimated mean squared errors, with each column corresponding to each of the regularization parameters in \code{lambda} and each row corresponding to each of the n.folds folds.}
#'  \item{mean.cv}{a nlambda-by-1 vector of the estimated mean squared errors, with each element corresponding to each of the regularization parameters in \code{lambda}.}
#'  \item{sd.cv}{a nlambda-by-1 vector of the standard deviation of the estimated mean squared errors, with each element corresponding to each of the regularization parameters in \code{lambda}.}
#'
#'
#' @author Park, Petkova, Tarpey, Ogden
#' @importFrom splines ns
#' @import SAM graphics stats
#' @seealso \code{samTEMsel}, \code{predict_samTEMsel}, \code{plot_samTEMsel}
#' @export
#'
#' @examples
#' set.seed(112)
#' n.train = 300
#' n.test  = 1000
#' n = n.train + n.test
#' p = 50
#' A = rbinom(n, 1, 0.5) + 1  # treatment variable taking a value in {1,2} with equal prob.
#' X = matrix(runif(n*p, -pi/2,pi/2), n, p)  # pretreatment covariates
#' noise = rnorm(n, 0, 0.5)
#' # X main effect on y; a highly nonlinear (cosine) function; depends on the first 10 covariates
#' main.effect = rep(0, n); for(j in 1:10){
#'   main.effect = main.effect + cos(X[,j])
#' }
#' # A-by-X ineraction effect on y; depends only on X1 and X2.
#' interaction.effect = (A-1.5)*X[,1]  +  2 * (A-1.5)*(cos(X[,2]) - 0.5)
#' # generate outcome y
#' y = main.effect  + interaction.effect + noise
#'
#'
#' # train/test set splitting
#' train.index = 1:n.train
#' y.train = y[train.index]
#' X.train = X[train.index,]
#' A.train = A[train.index]
#' y.test  = y[-train.index]
#' X.test  = X[-train.index,]
#' A.test  = A[-train.index]
#'
#'
#' # obtain an optimal regularization parameter by running cv.samTEMsel().
#' cv.obj = cv.samTEMsel(y.train, A.train, X.train, nlambda = 100)
#'
#' samTEMsel.obj = cv.obj$samTEMsel.obj
#' # samTEMsel.obj contains the sequence of fitted models over the grid of lambda
#' # see also, samTEMsel().
#'
#' # lambda.opt.index corresponds to the optimal regularization parameter chosen from cv.samTEMsel().
#' lambda.opt.index = cv.obj$lambda.opt.index
#' lambda.opt.index
#'
#' # plot the estimated component function of variable (j=)2, say.
#' plot_samTEMsel(samTEMsel.obj,  which.index = 2, lambda.index = lambda.opt.index)
#'
#' # make ITRs for subjects with pretreatment characteristics, X.test
#' trt.rule = make_ITR(samTEMsel.obj, newX = X.test, lambda.index = lambda.opt.index)$trt.rule
#' head(trt.rule)
#'
#' # an (IPWE) estimate of the "value" of this particualr treatment rule, trt.rule:
#' mean(y.test[A.test==trt.rule])
#'
#' # compare the above value to the following estimated "values" of "naive" treatment rules:
#' mean(y.test[A.test==1])   # a rule that assigns everyone to A=1
#' mean(y.test[A.test==2])   # a rule that assigns everyone to A=2
#'
cv.samTEMsel = function(y, A, X, mu.hat = NULL, d = 3, n.folds = 10, nlambda = 50, lambda.min.ratio = 0.01,
                        thol = 1e-05, max.ite = 1e+05, regfunc = "L1",
                        row.ordering = NULL, cv1sd = FALSE, terms.fit  = FALSE, plots=TRUE)
{

  # first fit the model based on the whole data, and obtain the fitted component functions
  samTEMsel.obj =  samTEMsel(y=y, A=A, X=X, d=d, terms.fit =terms.fit,
                             nlambda =nlambda, lambda.min.ratio=lambda.min.ratio,
                             thol=thol, max.ite=max.ite, regfunc=regfunc, mu.hat=mu.hat)

  # perform an n.folds cross-validation to select an optimal sparsity tuning parameter.
  if(is.null(row.ordering))  row.ordering = sample(length(y))   # randomly shuffle the row index
  datax   = data.frame(y=y, A=A, X=X)[row.ordering, ]
  basisx  = samTEMsel.obj$basis[row.ordering,  ]
  basiscx = samTEMsel.obj$basisc[row.ordering, ]
  folds   = cut(seq(1,nrow(datax)), breaks= n.folds, labels=FALSE)

  cv.storage =  matrix(NA, n.folds, nlambda)
  for(kk in 1:n.folds) {
    testIndexes =  which(folds==kk, arr.ind=TRUE)
    testData  = datax[testIndexes, ]
    trainData = datax[-testIndexes, ]
    train.fit = samTEMsel(y= trainData[,1], A =trainData[,2], X = trainData[,-c(1,2)], d=d, terms.fit =FALSE,
                          nlambda =nlambda, lambda.min.ratio=lambda.min.ratio,
                          thol=thol, max.ite=max.ite, regfunc=regfunc,
                          basis = basisx[-testIndexes,], basisc = basiscx[-testIndexes,], mu.hat= mu.hat[-testIndexes])
    test.pred =predict_samTEMsel(train.fit, newX=testData[,-c(1,2)], newA=testData[,2], basis = basisx[testIndexes,])

    for(l in 1:nlambda)  cv.storage[kk, l] = mean((testData[,1] - test.pred[,l])^2)
  }

  mean.cv = apply(cv.storage, 2, mean)
  sd.cv   = apply(cv.storage, 2, sd)/sqrt(n.folds)
  if(plots){
    plot(mean.cv, xlab=expression(paste("Regularization parameter ", lambda, " index")), ylab = "Mean squared error",cex.lab=1,lwd=2)
    lines(mean.cv+sd.cv, type = "l",  col ="red", lwd= 0.2)
    lines(mean.cv-sd.cv, type = "l",  col ="red", lwd= 0.2)
  }
  lambda.opt.index = which.min(mean.cv)
  if(cv1sd){
    cv1sd = mean.cv[lambda.opt.index] + sd.cv[lambda.opt.index] # move lambda.opt.index in the direction of increasing regularization
    lambda.opt.index = which(cv1sd > mean.cv)[1]
  }else{
    lambda.opt.index = which.min(mean.cv)
  }

  func_norm.opt = samTEMsel.obj$func_norm[, lambda.opt.index]
  nonzero.index = which(func_norm.opt > 0.1*max(func_norm.opt))

  return(list(samTEMsel.obj=samTEMsel.obj,
              lambda.opt.index=lambda.opt.index,
              nonzero.index=nonzero.index, func_norm.opt=func_norm.opt,
              cv.storage=cv.storage, mean.cv = mean.cv, sd.cv = sd.cv))
}





#' Sparse Additive Models for Treatment Effect-Modifier Selection (main function)
#'
#' The function \code{samTEMsel} implements estimation of a constrained sparse additve model.
#'
#' A constrained additive regression model represents the joint effects of treatment and p pretreatment covariates on an outcome via treatment-specific additive component functions defined over the p covariates, subject to the constraint that the expected value of the outcome given the covariates equals zero, while leaving the main effects of the covariates unspecified.
#' Under this flexible representation, the treatment-by-covariates interaction effects are determined by distinct shapes (across treatment levels) of the unspecified component functions.
#' Optimized under a penalized least square criterion with a L1 (or SCAD/MCP) penalty, the constrained additive model can effectively identify/select treatment effect-modifiers (from the p pretreatment covariates) that exhibit possibly nonlinear interactions with the treatment variable; this is achieved by producing a sparse set of estimated component functions.
#' The estimated nonzero component functions (available from the returned \code{samTEMsel} object) can be used to make individualized treatment recommendations (ITRs) for future subjects; see also \code{make_ITR} for such ITRs.
#'
#'
#' The regularization path is computed at a grid of values for the regularization parameter \code{lambda}.
#'
#'
#' @param y   a n-by-1 vector of responses
#' @param A   a n-by-1 vector of treatment variable; each element represents one of the L(>1) treatment conditions; e.g., c(1,2,1,1,3...); can be a factor-valued
#' @param X   a n-by-p matrix of pretreatment features
#' @param mu.hat   a n-by-1 vector of the fitted X main effect term of the model provided by the user; defult is \code{NULL}, in which case \code{mu.hat} is taken to be a vector of zeros; the optimal choice for this vector is E(y|X)
#' @param d     number of basis spline functions to be used for each component function; the default value is 3; d=1 corresponds to the linear model
#' @param thol  stopping precision; the default value is 1e-5.
#' @param max.ite  number of maximum iterations; the default value is 1e5.
#' @param regfunc  type of the regularizer; the default is "L1"; can also be "MCP" or "SCAD".
#' @param lambda   a user-supplied lambda sequence; typical usage is to have the program compute its own lambda sequence based on \code{nlambda} and \code{lambda.min.ratio}.
#' @param nlambda  total number of \code{lambda} values; the default value is 50.
#' @param lambda.min.ratio  the smallest value for lambda, as a fraction of lambda.max, the (data derived) entry value (i.e. the smallest value for which all coefficients are zero); the default is 0.01.
#' @param terms.fit   if \code{TRUE}, return the component-wise fitted vectors (\code{y.hat.list}) and the partial residuals (\code{resid.list}), evaluated over a grid of \code{lambda}.
#' @param basis   default is \code{NULL}; one can provide a n-by-d*p*L basis matrix associted with the spline representation of the model; this is to efficiently implement cross-validation in \code{cv.samTEMsel}.
#' @param basisc  default is \code{NULL}; one can provide a n-by-d*p*(L-1) basis matrix associted with the spline representation that incorporates the "orthogonality" constraint; this is only to efficiently implement cross-validation in \code{cv.samTEMsel}.
#'
#' @return a list of information of the fitted constrained additive models including
#'  \item{w}{the solution path matrix; the estimated d*p*L-by-nlambda coefficient matrix associted with B-splines, with each column corresponding to a regularization parameter in \code{lambda}}
#'  \item{wc}{the solution path matrix; the estimated d*p*(L-1)-by-nlambda coefficient matrix associted with B-splines that incorports the "orthogonality" constraint, with each column corresponding to a regularization parameter in \code{lambda}}
#'  \item{lambda}{a sequence of regularization parameter}
#'  \item{d}{the number of basis spline functions used for each component function}
#'  \item{func_norm}{the functional norm matrix (p-by-nlambda) with each column corresponding to a regularization parameter in \code{lambda}; since we have p variables, the length of each column is p.}
#'  \item{df}{the degree of freedom of the solution path (the number of non-zero component functions).}
#'  \item{knots}{a (d-1)-by-p matrix; each column contains the knots applied to the corresponding variable.}
#'  \item{Boundary.knots}{a 2-by-p matrix; each column contains the boundary points applied to the corresponding variable.}
#'  \item{sse}{sums of square errors of the fitted models over the solution path.}
#'  \item{intercept}{the list of L treatment-specific intercepts; each element of the list is a 1-by-nlambda matrix, with each column corresponding to a regularization parameter in \code{lambda}.}
#'  \item{X.min}{a p-by-1 vector, with each entry corresponding to the minimum of each input variable; used for rescaling in testing.}
#'  \item{X.ran}{a p-by-1 vector, with each entry corresponding to the range of each input variable; used for rescaling in testing.}
#'  \item{residuals}{the n-by-nlambda matrix of residuals of the fitted models, with each column corresponding to each regularization parameter in \code{lambda}.}
#'  \item{y.hat}{the n-by-nlambda matrix of fitted values for y, with each column corresponding to each regularization parameter in \code{lambda}.}
#'  \item{resid.list}{the list of p component-wise partial residuals, in which each (the jth) element is a n-by-nlambda matrix of the jth partial residuals, with each column corresponding to each regularization parameter in \code{lambda}.}
#'  \item{y.hat.list}{the list of p component-wise fitted values, in which each (the jth) element is a n-by-nlambda matrix of the jth fitted values, with each column corresponding to each regularization parameter in \code{lambda}.}
#'  \item{basis}{the n-by-d*p*L basis matrix associted with the spline representation of the (unconstrained) additive model.}
#'  \item{basisc}{the n-by-d*p*(L-1) basis matrix associted with the spline representation of the model that incorporates the "orthogonality" constraint.}
#'
#' @author Park, Petkova, Tarpey, Ogden
#' @import SAM
#' @seealso \code{cv.samTEMsel}, \code{predict_samTEMsel}, \code{plot_samTEMsel}, \code{make_ITR}
#' @export
#'
#' @examples
#'set.seed(112)
#' n.train = 400
#' n.test  = 50
#' n = n.train + n.test
#' p = 20
#' A = rbinom(n, 1, 0.5) + 1  # treatment variable taking a value in {1,2} with equal prob.
#' X = matrix(runif(n*p, -pi/2,pi/2), n, p)  # pretreatment covariates
#' noise = rnorm(n, 0, 0.5)
#' # X main effect on y; a highly nonlinear (cosine) function; depends on the first 10 covariates
#' main.effect = rep(0, n); for(j in 1:10){
#'   main.effect = main.effect + cos(X[,j])
#' }
#' # A-by-X ineraction effect on y; depends only on X1 and X2.
#' interaction.effect = (A-1.5)*X[,1]  +  2 * (A-1.5)*(cos(X[,2]) - 0.5)
#' # generate outcome y
#' y = main.effect  + interaction.effect + noise
#'
#'
#' # train/test set splitting
#' train.index = 1:n.train
#' y.train = y[train.index]
#' X.train = X[train.index,]
#' A.train = A[train.index]
#' y.test = y[-train.index]
#' X.test = X[-train.index,]
#' A.test = A[-train.index]
#'
#' # fit samTEMsel() based on the training set
#' samTEMsel.obj = samTEMsel(y.train, A.train, X.train, nlambda = 50)
#'
#' # a n.test-by-nlambda matrix of predicted values:
#' predict_samTEMsel(samTEMsel.obj, newX = X.test, newA = A.test)
#'
#' # pick a particular lambda.index, say, 10, as the regularization parameter.
#' lambda.index = 10
#' # for an optimal selection of lambda.index, see cv.samTEMsel().
#'
#' # a n.test-by-1 vector of predicted values (given lambda.index = 10)
#' predict_samTEMsel(samTEMsel.obj, newX = X.test, newA = A.test, lambda.index=lambda.index)
#'
#' # the estimated L2 norm of the component functions (given lambda.index=10)
#' samTEMsel.obj$func_norm[,lambda.index]
#'
#' # p component-wise fitted values (given lambda.index=10); the last column is the intercept
#' predict_samTEMsel(samTEMsel.obj, X.test, A.test, type="terms", lambda.index =lambda.index)
#'
#' # can plot the estimated component functions (say, the first two functions, j=1,2)
#' plot_samTEMsel(samTEMsel.obj, which.index = c(1,2), lambda.index = lambda.index)
#'
#' # can make inidividualized treatment recommendations (ITRs)
#' trt.rule = make_ITR(samTEMsel.obj, newX = X.test, lambda.index = lambda.index)$trt.rule
#' head(trt.rule)
#'
#' # an (IPWE) estimate of the "value" of this particualr treatment rule, trt.rule:
#' mean(y.test[A.test==trt.rule])
#'
#' # compare the above value to the following estimated "values" of "naive" treatment rules:
#' mean(y.test[A.test==1])   # just assign everyone to A=1
#' mean(y.test[A.test==2])   # just assign everyone to A=2
#'
samTEMsel =  function (y, A, X, mu.hat=NULL,
                       d = 3, # number of interior knots for b-spline expansion for each predictor
                       lambda = NULL, nlambda = 50, lambda.min.ratio = 0.01,
                       thol = 1e-05, max.ite = 1e+05, regfunc = "L1",
                       terms.fit = FALSE, basis= NULL, basisc= NULL)
{

  y = as.vector(y)
  A = as.numeric(as.factor(A))  # so that the variable A takes a value in 1, 2, 3,..
  X = as.matrix(X)
  # For discrete-valued variable, we jitter the values of the variables so that we can avaid sigular design matrices upon basis-expansion
  for(j in which(apply(X, 2, function(x) length(unique(x))) < d+2)){
    X[,j] = X[,j] + rnorm(length(X[,j]), 0, 0.01)
  }
  gcinfo(FALSE)
  fit = list()
  fit$d = d
  fit$X = X
  fit$A = A
  fit$y = y
  fit$mu.hat = mu.hat
  n = nrow(X)
  p = ncol(X)
  fit$p = p
  L = length(unique(A))
  A.unique = unique(A)
  fit$L = L
  fit$A.unique = A.unique

  A.indicator = matrix(0, n, L)
  for(a in 1:L)  A.indicator[,a] = as.numeric(A == A.unique[a])
  fit$A.indicator = A.indicator
  X.min = apply(X, 2, min)
  X.max = apply(X, 2, max)
  X.ran = X.max - X.min
  X.min.rep = matrix(rep(X.min, n), nrow = n, byrow = TRUE)
  X.ran.rep = matrix(rep(X.ran, n), nrow = n, byrow = TRUE)
  X.tilde   = (X - X.min.rep)/X.ran.rep
  fit$X.tilde = X.tilde
  fit$X.min = X.min
  fit$X.ran = X.ran

  knots = matrix(0, d - 1, p)
  Boundary.knots = matrix(0, 2, p)

  compute.basis = 0
  if(is.null(basis)){
    compute.basis  = 1
    basis  = matrix(0, n, d*p*L)       # basis will be the design matrix
  }
  if(is.null(basisc)) basisc = matrix(0, n, d*p*(L-1))   # basisc will absorb the linear constraint into its construction.

  pr.A = summary(as.factor(A))/n   # probability of the treatment assignment
  C = NULL  # the constraint on the coefficients to be encoded as a matrix
  for(a in 1:L)  C = cbind(C, diag(rep(pr.A[a], d)))
  if(d==1) C = matrix(pr.A, 1, L)
  qrc = qr(t(C))
  null.space.C = qr.Q(qrc,complete=TRUE)[,(nrow(C)+1):ncol(C)]

  if(compute.basis){
    for (j in 1:p){
      index = (j - 1) * (d*L) + c(1:(d*L))
      indexc = (j - 1) * (d*(L-1)) + c(1:(d*(L-1)))
      tmp0 = ns(X.tilde[, j], df = d)
      knots[, j] = attr(tmp0, "knots")
      Boundary.knots[, j] = attr(tmp0, "Boundary.knots")
      tmp1 = NULL
      for(a in 1:L)  tmp1 = cbind(tmp1, A.indicator[,a]*tmp0)
      basis[, index]   = tmp1
      basisc[, indexc] = tmp1 %*% null.space.C
    }
  }

  basisc.mean = apply(basisc, 2, mean)
  basisc.mean.rep = matrix(rep(basisc.mean, n), nrow = n, byrow = TRUE)
  basisc = basisc - basisc.mean.rep
  fit$basis  = basis
  fit$basisc = basisc
  fit$knots  = knots
  fit$Boundary.knots = Boundary.knots

  y.mean  = vector(mode = "list", length = L)
  if(is.null(mu.hat)){y.tilde = y} else y.tilde = y - mu.hat
  for(a in 1:L){
    y.mean[[a]] = mean(y.tilde[A==A.unique[a]])
    y.tilde = y.tilde - y.mean[[a]]*(A==A.unique[a])
  }
  fit$y.tilde = y.tilde  # now, y.tilde does not have the X and A main effects in it

  lambda_input = 1
  if (is.null(lambda)){
    lambda_input = 0
    if (is.null(nlambda))
      nlambda = 30
    lambda = exp(seq(log(1), log(lambda.min.ratio), length = nlambda))
  }else nlambda = length(lambda)

  out = .C("grplasso", y = as.double(y.tilde), X = as.double(basisc),
           lambda = as.double(lambda), nnlambda = as.integer(nlambda),
           nn = as.integer(n), dd = as.integer(p), pp = as.integer(d*(L-1)),
           ww = as.double(matrix(0, d * p*(L-1), nlambda)), mmax_ite = as.integer(max.ite),
           tthol = as.double(thol), regfunc = as.character(regfunc),
           iinput = as.integer(lambda_input), df = as.integer(rep(0, nlambda)),
           sse = as.double(rep(0, nlambda)), func_norm = as.double(matrix(0, p, nlambda)),
           PACKAGE = "SAM")
  fit$lambda = out$lambda
  fit$wc = matrix(out$w, ncol = nlambda)
  fit$w  = matrix(0, d * p *L, nlambda)
  for (j in 1:p){
    index = (j - 1) * (d*L) + c(1:(d*L))
    indexc = (j - 1) * (d*(L-1)) + c(1:(d*(L-1)))
    fit$w[index,] = as.matrix(null.space.C) %*% fit$wc[indexc,]
  }

  fit$df = out$df
  fit$sse = out$sse
  fit$func_norm = matrix(out$func_norm, ncol = nlambda)
  fit$intercept = vector(mode = "list", length = L)
  for(a in 1:L) fit$intercept[[a]] = rep(y.mean[[a]], nlambda) - t(basisc.mean) %*% fit$wc
  intercepts = NULL
  for(a in 1:L) intercepts = rbind(intercepts, fit$intercept[[a]])
  fit$y.hat = cbind(basis, A.indicator) %*% rbind(as.matrix(fit$w), intercepts)

  # Compute partial residuals
  y.hat.list  = resid.list =  vector(mode = "list", length = p)
  fit$residuals  = y.tilde - basis %*% as.matrix(fit$w)
  if(terms.fit){
    for (j in 1:p){
      index = (j - 1) * (d*L) + c(1:(d*L))
      y.hat.list[[j]] = basis[,index] %*% as.matrix(fit$w[index, ])
      resid.list[[j]] = fit$residuals + y.hat.list[[j]]
    }
    fit$y.hat.list = y.hat.list
    fit$resid.list = resid.list
  }

  rm(out, X, X.tilde, y, y.tilde, basis, basisc, X.min.rep, X.ran.rep, basisc.mean.rep)

  class(fit) = "samTEMsel"
  return(fit)
}





#' \code{samTEMsel} prediction function
#'
#' \code{predict_samTEMsel} makes predictions given a (new) set of covariates \code{newX} and a (new) vector of treatment indicators \code{newA} based on a constrained sparse additive model \code{samTEMsel.obj}.
#' Specifically, \code{predict_samTEMsel} predicts the responses y based on the X-by-A interaction effect (and the A main effect) portion of the model.
#'
#' @param samTEMsel.obj  a \code{samTEMsel} object
#' @param newX a (n by p) matrix of new values for the covariates X at which predictions are to be made; if \code{NULL}, X from the training set is used.
#' @param newA a (n by 1) vector of new values for the treatment A at which predictions are to be made; if \code{NULL}, A from the training set is used.
#' @param type the type of prediction required; the default "response" gives the predicted responses y based on the whole model; the alternative "terms" gives the component-wise predicted responses from each of the p components (and plus the treatment-specific intercepts) of the model.
#' @param lambda.index an index of the tuning parameter \code{lambda} at which predictions are to be made; one can supply \code{lambda.opt.index} obtained from the function \code{cv.samTEMsel}; the default is \code{NULL}, in which case the predictions based on the most non-sparse model is returned.
#' @param basis  a basis (design) matrix associated with the testing set (newX and newA) provided by the user; the default is \code{NULL}; this is only to efficiently implement the cross-validation in \code{cv.samTEMsel}.
#'
#' @return
#' \item{value}{a (n-by-\code{length(lambda.index)}) matrix of predicted values; a (n-by-\code{length(lambda.index)}*(p+1)) matrix of predicted values if type =  "terms".}
#'
#'
#' @author Park, Petkova, Tarpey, Ogden
#' @seealso \code{samTEMsel},\code{cv.samTEMsel}, \code{plot_samTEMsel}
#' @export
#'
#'
predict_samTEMsel = function(samTEMsel.obj, newX=NULL, newA=NULL,
                             type = "response", # "terms"
                             lambda.index = NULL,
                             basis= NULL)
{

  if(!inherits(samTEMsel.obj, "samTEMsel"))   # checks input
    stop("samTEMsel.obj must be of class `samTEMsel'")

  nlambda = length(samTEMsel.obj$lambda)
  if (is.null(lambda.index)) {
    lambda.index = seq(nlambda)
  }else {
    lambda.index = intersect(lambda.index, seq(nlambda))
  }

  gcinfo(FALSE)
  out = list()
  if(is.null(newX)) newX = samTEMsel.obj$X
  if(is.null(newA)) newA = samTEMsel.obj$A
  n = nrow(as.matrix(newX))
  X.min.rep = matrix(rep(samTEMsel.obj$X.min, n), nrow=n, byrow=TRUE)
  X.ran.rep = matrix(rep(samTEMsel.obj$X.ran, n), nrow=n, byrow=TRUE)
  X.tilde   = (as.matrix(newX) - X.min.rep)/X.ran.rep
  X.tilde[X.tilde < 0] = 0.01
  X.tilde[X.tilde > 1] = 0.99

  L = samTEMsel.obj$L
  A.unique    = unique(samTEMsel.obj$A)
  A.indicator = matrix(0, n, L)
  for(a in 1:L) A.indicator[,a] = as.numeric(newA == A.unique[a])
  d = samTEMsel.obj$d
  p = samTEMsel.obj$p
  compute.basis = 0
  if(is.null(basis)){
    compute.basis  = 1
    basis  = matrix(0, n, d*p*L)       # basis will be the model matrix
  }

  if(compute.basis){
    for(j in 1:p){    # construct the jth model matrix
      index = (j - 1) * d*L + c(1:(d*L))
      tmp0 = ns(X.tilde[,j], df= d, knots= samTEMsel.obj$knots[,j], Boundary.knots= samTEMsel.obj$Boundary.knots[,j])
      tmp1 = NULL
      for(a in 1:L){
        tmp1 = cbind(tmp1, tmp0*A.indicator[,a])
      }
      basis[,index] = tmp1
    }
  }

  # compute the treatment-specific intercepts
  intercepts =  NULL
  for(a in 1:L){
    intercepts =  rbind(intercepts, samTEMsel.obj$intercept[[a]][, lambda.index])
  }

  # make predictions
  if(type=="terms"){
    y.hat.list = vector(mode = "list", length = p)
    for (j in 1:p){
      index = (j - 1) *d*L + c(1:(d*L))
      y.hat.list[[j]] = basis[,index] %*% as.matrix(samTEMsel.obj$w[index, lambda.index])
    }
    values = cbind(do.call("cbind", y.hat.list), A.indicator %*% intercepts)
  }else{
    values = cbind(basis, A.indicator) %*% rbind(as.matrix(samTEMsel.obj$w[, lambda.index]), intercepts)
  }

  return(values)
}




#' make individualized treatment recommendations (ITRs) based on a \code{samTEMsel} object
#'
#' The function \code{make_ITR} returns individualized treatment decision recommendations for subjects with pretreatment characteristics \code{newX}, given a \code{samTEMsel} object, \code{samTEMsel.obj}, and an (optimal) regularization parameter index, \code{lambda.index}.
#'
#'
#' @param samTEMsel.obj  a \code{samTEMsel} object, containing the fitted models.
#' @param newX a (n-by-p) matrix of new values for the covariates X at which predictions are to be made; if \code{NULL}, X from the training set is used.
#' @param lambda.index an index of the regularization parameter \code{lambda} at which predictions are to be made; one can supply \code{lambda.opt.index} obtained from the function \code{cv.samTEMsel()}; the default is \code{NULL}, in which case the predictions are made based on the most non-sparse model.
#' @param maximize default is \code{TRUE}, assuming a larger value of the outcome is better; if \code{FALSE}, a smaller value is assumed to be prefered.
#'
#' @return
#' \item{pred.new}{a (n-by-L) matrix of predicted values, with each column representing one of the L treatment options.}
#' \item{trt.rule}{a (n-by-1) vector of the individualized treatment recommendations}
#'
#' @author Park, Petkova, Tarpey, Ogden
#' @seealso \code{samTEMsel},\code{cv.samTEMsel}, \code{predict_samTEMsel}
#' @export
#'
# make treatment decisions given newX, based on the estimated sequence of models, samTEMsel.obj, and a given regularization parameter index, lambda.index.
make_ITR  = function(samTEMsel.obj, newX = NULL, lambda.index = NULL, maximize=TRUE){

  if(!inherits(samTEMsel.obj, "samTEMsel"))   # checks input
    stop("samTEMsel.obj must be of class `samTEMsel'")

  if (is.null(lambda.index))  lambda.index = length(samTEMsel.obj$lambda)

  L = samTEMsel.obj$L
  newX = as.matrix(newX)
  n = nrow(newX)
  # compute potential outcome for each of the L treatment conditions, given newX.
  pred.new= matrix(0, nrow=n, ncol = L)
  for(a in 1:L){
    pred.new[,a] =predict_samTEMsel(samTEMsel.obj, newX= newX,  newA=rep(a, n), lambda.index = lambda.index)
  }
  # compute treatment assignment
  if(maximize){
    trt.rule = apply(pred.new, 1, which.max)
  }else{
    trt.rule = apply(pred.new, 1, which.min)
  }
  if(L==2)  colnames(pred.new) <- c("Tr1", "Tr2")

  return(list(trt.rule = trt.rule, pred.new = pred.new))
}





#' plot component functions from a \code{samTEMsel} object
#'
#' Produces plots of the component functions from a \code{samTEMsel} object.
#'
#' @param samTEMsel.obj  a \code{samTEMsel} object
#' @param newX a (n-by-p) matrix of new values for the covariates X at which plots are to be made; the default is \code{NULL}, in which case X is taken from the training set.
#' @param newA a (n-by-1) vector of new values for the treatment A at which plots are to be made; the default is \code{NULL}, in which case A is taken from the training set.
#' @param scatter.plot if \code{TRUE}, draw scatter plots of partial residuals versus the covariates; these scatter plots are made based on the training observations;  the default is \code{TRUE}.
#' @param lambda.index an index of the tuning parameter \code{lambda} at which plots are to be made; one can supply \code{lambda.opt.index} obtained from the function \code{cv.samTEMsel}; the default is \code{NULL}, in which case \code{plot_samTEMsel} utilizes the most non-sparse model.
#' @param which.index  this specifies which component functions are to be plotted; the default is all p component functions, i.e., 1:p.
#' @param ylims this specifies the vertical range of the plots, e.g., c(-10, 10).
#' @param solution.path  if \code{TRUE}, draw the functional norms of the fitted component functions (based on the training set) versus the regularization parameter; the default is \code{FALSE}.
#'
#'
#' @author Park, Petkova, Tarpey, Ogden
#' @seealso \code{samTEMsel},\code{predict_samTEMsel}, \code{cv.samTEMsel}
#' @export
#'
#'
plot_samTEMsel <- function(samTEMsel.obj, newX=NULL, newA=NULL, scatter.plot = TRUE,
                           lambda.index, which.index, ylims, solution.path = FALSE)
{

  if(!inherits(samTEMsel.obj, "samTEMsel"))   # checks input
      stop("samTEMsel.obj must be of class `samTEMsel'")
  p = samTEMsel.obj$p
  L = samTEMsel.obj$L

  if(missing(which.index))  which.index = 1:p
  if(missing(lambda.index)) lambda.index = length(samTEMsel.obj$lambda)

  # compute component-wise predictions, given newX and newA
  predicted =predict_samTEMsel(samTEMsel.obj, newX=newX, newA=newA, type = "terms", lambda.index = lambda.index)

  if(scatter.plot){
    resid.list =  vector(mode = "list", length = p)
    predicted.train = predict_samTEMsel(samTEMsel.obj, type = "terms", lambda.index = lambda.index)
    for(j in 1:p){
      resid.list[[j]] = samTEMsel.obj$residuals[,lambda.index] + predicted.train[,j]
    }
  }
  if(missing(ylims)){
    if(scatter.plot){
      ylims =  range(do.call("cbind", resid.list))
    }else{
      ylims = range(predicted[, which.index])
    }
  }
  if(is.null(newX)) newX = samTEMsel.obj$X
  if(is.null(newA)) newA = samTEMsel.obj$A

  colvals    =  1:L
  colvals[1] = "royalblue3"
  colvals[L] = "tomato3"
  if(L > 2)  colvals[2] = "darkseagreen"

  for(j in which.index){
    o = order(newX[newA==1, j])
    plot(matrix(newX[newA==1,],ncol=p)[o, j], predicted[newA==1, j][o], type = "l",
         ylab = paste("f(X", j, ")", sep = ""), xlab = paste("X", j, sep = ""),
         ylim = ylims, xlim= c(range(newX[,j])[1], range(newX[,j])[2]),
         col = colvals[1], lwd = 2)
    if(scatter.plot)  points(samTEMsel.obj$X[samTEMsel.obj$A==1,j],
                             resid.list[[j]][samTEMsel.obj$A==1],
                             col = colvals[1])
    if(L > 1){
      for(a in 2:L){
        o = order(newX[newA==a, j])
        lines(matrix(newX[newA==a,],ncol=p)[o, j], predicted[newA==a, j][o], col = colvals[a], lwd = 2)
        if(scatter.plot)  points(samTEMsel.obj$X[samTEMsel.obj$A==a,j],
                                 resid.list[[j]][samTEMsel.obj$A==a],
                                 col = colvals[a], pch = a, lwd = 1)
      }
    }
  }

  if(solution.path){
    par = par(omi = c(0.0, 0.0, 0, 0), mai = c(1, 1, 0.1, 0.1))
    matplot(samTEMsel.obj$lambda[length(samTEMsel.obj$lambda):1],t(samTEMsel.obj$func_norm),type="l",
            xlab=expression(paste("Regularization parameter ", lambda)),ylab = "Functional norms",cex.lab=1,log="x",lwd=2)
  }

}



######################################################################
## END OF THE FILE
######################################################################
