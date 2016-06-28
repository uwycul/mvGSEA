edited_speedglm_L1 <- function (formula, data, family = gaussian(), weights = NULL, 
                                start = NULL, etastart = NULL, mustart = NULL, offset = NULL, 
                                maxit = 25, k = 2, sparse = NULL, set.default = list(), trace = FALSE, 
                                method = c("eigen", "Cholesky", "qr"), model = FALSE, y = FALSE, 
                                fitted = FALSE, ...) 
{
  data <- environment(formula)
  family =
  call <- match.call()
  target <- y
  M <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset"), names(M), 0L)
  M <- M[c(1L, m)]
  M$drop.unused.levels <- TRUE
  M[[1L]] <- quote(stats::model.frame)
  M <- eval(M, parent.frame())
  y <- M[[1]]
  tf <- terms(formula, data = data)
  X <- model.matrix(tf, M)
  offset <- model.offset(M)
  intercept <- attributes(tf)$intercept
  set <- list(sparselim = 0.9, camp = 0.01, eigendec = TRUE, 
              row.chunk = NULL, tol.solve = .Machine$double.eps, acc = 1e-08, 
              tol.values = 1e-07, tol.vectors = 1e-07, method = match.arg(method))
  nmsC <- names(set)
  set[(namc <- names(set.default))] <- set.default
  if (length(noNms <- namc[!namc %in% nmsC]) > 0) 
    warning("unknown names in set.default: ", paste(noNms, 
                                                    collapse = ", "))
  rval <- edited_L1_speedglm.wfit(y = y, X = X, family = family, weights = weights, 
                        start = start, etastart = etastart, mustart = mustart, 
                        offset = offset, intercept = intercept, row.chunk = set$row.chunk, 
                        maxit = maxit, k = k, acc = set$acc, sparselim = set$sparselim, 
                        camp = set$camp, eigendec = set$eigendec, tol.solve = set$tol.solve, 
                        sparse = sparse, tol.values = set$tol.values, trace = trace, 
                        tol.vectors = set$tol.vectors, method = set$method)
  rval$terms <- tf
  rval$call <- call
  class(rval) <- c("speedglm", "speedlm")
  if (model) 
    rval$model <- M
  if (fitted) 
    rval$linear.predictors <- predict.speedlm(rval, newdata = M)
  if (target) 
    rval$y <- y
  if ((rval$iter == maxit) & (!rval$convergence)) 
    warning("Maximum number of iterations reached without convergence")
  rval
}