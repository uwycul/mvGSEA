
#contains collection of functions necessary for the speedglm computation

ll.speedglm <- function(family,aic.model,nvar){
  switch(family,
         binomial = -(aic.model-2 * nvar)/2)
}

predict.speedglm <- function (object, newdata, type = c("link", "response"), na.action = na.pass, 
          ...) 
{
  family=binomial()
  type <- match.arg(type)
  if (missing(newdata) & is.null(object$linear.predictors)) 
    warning("fitted values were not returned from the speedglm object: \n            use the original data by setting argument 'newdata' or refit \n            the model by specifying fitted=TRUE.")
  na.act <- object$na.action
  object$na.action <- NULL
  if (missing(newdata)) {
    pred <- switch(type, link = object$linear.predictors, 
                   response = fitted(object))
    if (!is.null(na.act)) 
      pred <- napredict(na.act, pred)
  }
  else {
    pred <- predict.speedlm(object, newdata, type = "response", 
                            na.action = na.action)
    switch(type, response = {
      pred <- family$linkinv(pred)
    }, link = )
  }
  pred
}

predict.speedlm <- function (object, newdata, na.action = na.pass, ...) 
{
  family=binomial()
  tt <- terms(object)
  if (!inherits(object, c("speedlm", "edited_speedglm_L1"))) 
    warning("calling predict.speedlm(<fake-speedlm/speedglm-object>) ...")
  if (missing(newdata) || is.null(newdata)) {
    if (is.null(object$fitted.values)) 
      warning("fitted values were not returned from the speedglm object: \n              use the original data by setting argument 'newdata' or refit \n              the model by specifying fitted=TRUE.")
    return(object$fitted.values)
  }
  else {
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.action)
    if (!is.null(cl <- attr(Terms, "dataClasses"))) 
      .checkMFClasses(cl, m)
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    offset <- rep(0, nrow(X))
    if (!is.null(off.num <- attr(tt, "offset"))) 
      for (i in off.num) offset <- offset + eval(attr(tt, 
                                                      "variables")[[i + 1]], newdata)
    if (!is.null(object$call$offset)) 
      offset <- offset + eval(object$call$offset, newdata)
  }
  p <- object$rank
  ord <- colnames(X)
  if (p < ncol(X) && !(missing(newdata) || is.null(newdata))) 
    warning("prediction from a rank-deficient fit may be misleading")
  beta <- object$coefficients
  beta[is.na(beta)] <- 0
  predictor <- drop(X[, ord, drop = FALSE] %*% beta[ord])
  if (!is.null(offset)) 
    predictor <- predictor + offset
  if (missing(newdata) && !is.null(na.act <- object$na.action)) 
    predictor <- napredict(na.act, predictor)
  predictor
}

coef.speedglm <- function (object, ...) {
  object$coefficients
}

vcov.speedglm<-{
function (object, ...) 
  object$dispersion * solve(object$XTX)
  }



aic.binomial<-function(y, n, mu, wt, dev){
  m <- if (any(n > 1)) n else wt
  -2 * sum(ifelse(m > 0, (wt/m), 0) * dbinom(round(m * y),
                                             round(m), mu, log = TRUE))
}

is.sparse <- function (X, sparselim = 0.9, camp = 0.05) 
{
  if (prod(dim(X)) < 500) 
    camp <- 1
  subX <- sample(X, round((nrow(X) * ncol(X) * camp), digits = 0), 
                 replace = FALSE)
  p <- sum(subX == 0)/prod(length(subX))
  if (p > sparselim) 
    sparse <- TRUE
  else sparse <- FALSE
  attr(sparse, "proportion of zeros in the sample") <- p
  sparse
}

cp <- function (X, w = NULL, row.chunk = NULL, sparse = FALSE) 
{
  if (sparse) {
    if (!is.null(w)) 
      X <- sqrt(w) * X
    if (class(X) != "dgCMatrix") 
      X <- as(X, "dgCMatrix")
    new.B <- crossprod(X)
  }
  else {
    if (is.null(row.chunk)) {
      new.B <- if (is.null(w)) 
        crossprod(X)
      else crossprod(sqrt(w) * X)
    }
    else {
      if (row.chunk >= (nrow(X) - 1)) 
        row.chunk <- nrow(X)
      if (row.chunk <= 1) 
        row.chunk <- 2
      mod <- nrow(X)%%row.chunk
      last.block <- (mod > 0)
      G <- nrow(X)%/%row.chunk - (last.block * (2 * row.chunk <= 
                                                  nrow(X)))
      new.B <- matrix(0, ncol(X), ncol(X))
      a <- row.chunk
      j <- 1
      if (is.null(w)) {
        for (i in 1:G) {
          B <- crossprod(X[j:(i * a), ])
          new.B <- new.B + B
          j <- j + a
        }
        if (mod > 0) 
          new.B <- new.B + crossprod(X[j:nrow(X), ])
      }
      else {
        for (i in 1:G) {
          B <- crossprod(sqrt(w[j:(i * a)]) * X[j:(i * 
                                                     a), ])
          new.B <- new.B + B
          j <- j + a
        }
        if (mod > 0) 
          new.B <- new.B + crossprod(sqrt(w[j:nrow(X)]) * 
                                       X[j:nrow(X), ])
      }
    }
  }
  as(new.B, "matrix")
}

control <- function (B, symmetric = TRUE, tol.values = 1e-07, tol.vectors = 1e-07, 
                     out.B = TRUE, method = c("eigen", "Cholesky")) 
{
  method <- match.arg(method)
  if (!(method %in% c("eigen", "Cholesky"))) 
    stop("method not valid or not implemented")
  if (method == "eigen") {
    n <- ncol(B)
    sa <- 1:n
    nok <- NULL
    auto <- eigen(B, symmetric, only.values = TRUE)
    totcoll <- sum(abs(auto$values) < tol.values)
    ncoll <- totcoll
    rank <- n - ncoll
    i <- 1
    while (ncoll != 0) {
      auto <- eigen(B, symmetric)
      j <- as.matrix(abs(auto$vectors[, n]) < tol.vectors)
      coll <- which(!j)
      coll <- coll[length(coll)]
      B <- B[-coll, -coll]
      nok[i] <- coll
      ncoll <- sum(abs(auto$values) < tol.values) - 1
      n <- ncol(B)
      i <- i + 1
    }
    ok <- if (!is.null(nok)) 
      sa[-nok]
    else sa
  }
  if (method == "Cholesky") {
    A <- chol(B, pivot = TRUE)
    pivot <- attributes(A)$pivot
    rank <- attributes(A)$rank
    ok <- sort(pivot[1:rank])
    nok <- if (rank < length(pivot)) 
      pivot[(rank + 1):length(pivot)]
    else NULL
    B <- B[ok, ok]
  }
  rval <- if (out.B) 
    list(XTX = B, rank = rank, pivot = c(ok, nok))
  else list(rank = rank, pivot = c(ok, nok))
  rval
}
