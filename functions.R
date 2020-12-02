
options(warn=-1, digits=4)

orthoBasis = function(order, denseGrid, type, normalized = T){
  if (!("orthopolynom" %in% rownames(installed.packages()))) 
    install.packages("orthopolynom")
  library(orthopolynom)
  
  if (type == 'shiftedLegendre') {
    poly.value = polynomial.values(slegendre.polynomials(n = max(order), normalized = normalized), x = denseGrid)
    res = matrix(NA, nrow = length(order), ncol = length(denseGrid))
    for (i in 1:length(order)){
      res[i, ] = poly.value[[order[i]+1]]
    }
  }
  return(res)
}

crossVali = function(y.old, x.old, y.new, x.new, timePts, nfold = 5, method = c('PCC', 'PLCC'), pars){

  parsOpt = rep(list(NULL), 2)
  names(parsOpt) = c('PCC', 'PLCC')
  
  locOpt = rep(list(NULL), 2)
  names(parsOpt) = c('PCC', 'PLCC')

  cvVec = rep(list(NULL), 2)
  names(cvVec) = c('PCC', 'PLCC')
  
  errPred = rep(list(NULL), 2)
  names(cvVec) = c('PCC', 'PLCC')

  n = nrow(x.old)
  test = split(1:n, 1:nfold)

  for (k in 1:nfold) {
    x.train = x.old[-test[[k]],]
    y.train = y.old[-test[[k]]]
    x.test = x.old[test[[k]],]
    y.test = y.old[test[[k]]]

    if ('PCC' %in% method) {
      resultPCC = PCC(y.train, x.train, y.test, x.test, timePts, pars$PCC)
      if (k == 1) cvVec$PCC = resultPCC$errPredVec
      else cvVec$PCC = resultPCC$errPredVec + cvVec$PCC
    }
    if ('PLCC' %in% method) {
      resultPLCC = PLCC(y.train, x.train, y.test, x.test, timePts, pars$PLCC)
      if (k == 1) cvVec$PLCC = resultPLCC$errPredVec
      else cvVec$PLCC = resultPLCC$errPredVec + cvVec$PLCC
    }
  }
  
  if (any(method %in% c('PCC'))) {
    locOpt$PCC = which(cvVec$PCC == min(cvVec$PCC, na.rm = TRUE))
    parsOpt$PCC = matrix(pars$PCC[ , locOpt$PCC])
    resultPCC = PCC(y.old, x.old, y.new, x.new, timePts, parsOpt$PCC)
    errPred$PCC = resultPCC$errPredVec
  }

  if (any(method %in% c('PLCC'))) {
    locOpt$PLCC = which(cvVec$PLCC == min(cvVec$PLCC, na.rm = TRUE))
    parsOpt$PLCC = matrix(pars$PLCC[locOpt$PLCC])
    resultPLCC = PLCC(y.old, x.old, y.new, x.new, timePts, parsOpt$PLCC)
    errPred$PLCC = resultPLCC$errPredVec
  }
  return(list(parsOpt = parsOpt, errPred = errPred))
}

#' compute pMax
#'
#'

pUpper.compu = function(x, y, timePts, proportion){
  if (!("fda" %in% rownames(installed.packages()))) install.packages("fda")
  library(fda)

  index.0 = (y==0)
  index.1 = (y==1)

  mu.0 = colMeans(x[index.0,])
  mu.1 = colMeans(x[index.1,])

  expansion.x = expand.bspline(x, timePts = timePts)
  W = inprod(expansion.x$basis, expansion.x$basis)
  x.fd = smooth.basisPar(argvals=timePts,
                         y = t(rbind(sweep(x[index.0,], 2, mu.0), sweep(x[index.1,], 2, mu.1))),
                         fdobj=expansion.x$basis,
                         Lfdobj=int2Lfd(2),
                         lambda=expansion.x$lambda)
  x.pca.obj = pca.fd(x.fd$fd, nharm=min(dim(x)-1), centerfns=F)
  pMax = which.max(cumsum(x.pca.obj$varprop)>=proportion)
  return(pMax)
}

#'
#' 
#'   

impute.curve = function(x){
  if (!("spatialEco" %in% rownames(installed.packages()))) install.packages("spatialEco")
  library(spatialEco)
  
  gcv = NULL
  lambda.vec = 10^(0:6)
  
  for (lambda in lambda.vec){
    x.smooth = array(NA, dim = dim(x))
    for (i in 1:nrow(x)){
      x.smooth[i, ] = poly.regression(x[i, ], s = lambda)
    }
    gcv = c(gcv, sum(((x-x.smooth)[!is.na(x-x.smooth)])^2))
  }
  lambda = lambda.vec[which.min(gcv)]
  for (i in 1:nrow(x)){
    x.smooth[i, ] = poly.regression(x[i, ], s = lambda, impute = T, na.only = T)
  }
  
  return(x.smooth)
}

#' expand.bspline is to expand 1-dim functional data with respect to a B-spline basis.
#' The output martrix C is guaranteed to have full column rank.
#' @param x is a matrix. Each row represents a sample.
#' @param nOrder=4 by default indicates the cubic B-spline.

expand.bspline = function(x, timePts = NULL, nOrder = 4){
  if (!("fda" %in% rownames(installed.packages()))) install.packages("fda")
  if (!all(c(
    exists("create.bspline.basis", mode="function"),
    exists("int2Lfd", mode="function"),
    exists("fdPar", mode="function"),
    exists("Data2fd", mode="function"),
    exists("inprod", mode="function"),
    exists('deriv.fd', mode='function')
  ))) library("fda")
  
  if (is.null(timePts)) timePts = 1:ncol(x)
  K = min(nOrder+ncol(x)-2, nrow(x)-1) # number of basis functions
  splineBasis = create.bspline.basis(rangeval=range(timePts), nbasis=K, norder=nOrder, names="bspl")
  D2Lfd = int2Lfd(m=2)
  
  # tune the value of lambda for smoothing x
  log10lambda = 0:6
  gcvsave = NULL
  for (i in log10lambda) {
    lambda = 10^i
    D2fdPar = fdPar(splineBasis, D2Lfd, lambda)
    smooth = smooth.basis(y=t(x), argvals=timePts, D2fdPar, dfscale=1.2)
    gcvsave = c(gcvsave, sum(smooth$gcv[!is.nan(smooth$gcv)]))
  }
  
  lambda = 10^log10lambda[which.min(gcvsave)]
  D2fdPar = fdPar(splineBasis, D2Lfd, lambda)
  smooth = smooth.basis(y=t(x), argvals=timePts, D2fdPar, dfscale=1.2)
  
  return(list(fdObj = smooth$fd, coef = t(smooth$fd$coefs), basis = splineBasis, D2fdPar = D2fdPar, lambda = lambda))
} 

#' try to find the constraint global minimizer of f
#' 
#' @param ninitial is the number of initial values.

manymodes = function(f, ninitial){
  
  initial1 = -1/(ninitial+1)*(1:ninitial)
  initial2 = 1e8/(ninitial+1)*(1:ninitial)
  
  for (i in 1:length(initial1)) {
    result.tmp = optim(initial1[i], f, method=c("L-BFGS-B"), lower=-1+epsilon, upper=-epsilon)
    if (i==1) result = result.tmp
    if (result.tmp$value<result$value) result = result.tmp
  }
  for (i in 1:length(initial2)) {
    result.tmp = optim(initial2[i], f, method=c("L-BFGS-B"), lower=epsilon, upper=Inf)
    if (result.tmp$value<result$value) result = result.tmp
  }
  
  x.star = result$par
  return(x.star)
}

#'
#'
#'

half = function(A){
  svdA = svd(A)
  E.half = diag(sqrt(svdA$d))
  U = svdA$u
  V = svdA$v
  return(U %*% E.half %*% t(V))
}

#'
#'
#'

halfinv = function(A){
  svdA = svd(A)
  tmp = svdA$d[abs(svdA$d) > epsilon]
  E.halfinv = diag(c(1/sqrt(tmp), rep(0, times=length(svdA$d)-length(tmp))))
  U = svdA$u
  V = svdA$v
  return(U %*% E.halfinv %*% t(V))
}

#' SIMPLS by DeJong(1993) pp.262
#' @input n*p matrix x, n*m matrix y, number of components
#'

simpls = function(y, x, ncomp){
  
  if(is.vector(y)) y = matrix(y) # transform vector y into a column vector
  
  Rmat = NULL
  Tmat = NULL
  Pmat = NULL
  Qmat = NULL
  Umat = NULL
  Vmat = NULL
  
  if (ncol(y) > 1) Y0 = sweep(y, 2, colMeans(y))
  else Y0 = y - mean(y)
  X = x
  
  S = crossprod(X, Y0)
  
  for (a in 1:ncomp){
    q = svd(crossprod(S), nu=1, nv=1)$u # Y block factor weights
    r = S %*% q # X block factor weights
    t = X %*% r # X block factor scores
    t = t - mean(t) # center scores
    normt = sum(t^2)^.5
    t = t/normt # normalize scores
    r = r/normt # adapt weights a~~ordingIy
    p = crossprod(X, t) # X block factor loadings
    q = crossprod(Y0, t) # Y block factor loadings
    u = Y0 %*% q # Y block factor scores
    v = p # initialize orthogonal loadings
    if (a>1){
      v = v - Vmat %*% crossprod(Vmat, p) # make v perpendicular of previous loadings
      u = u - Tmat %*% crossprod(Tmat, u) # make u perpendicular of previous t values
    }
    v = v / sum(v^2)^.5 # normalize orthogonal loadings
    S = S - v %*% crossprod(v, S) # deflate S with respect to current loadings
    
    Rmat = cbind(Rmat, r)
    Tmat = cbind(Tmat, t)
    Pmat = cbind(Pmat, p)
    Qmat = cbind(Qmat, q)
    Umat = cbind(Umat, u)
    Vmat = cbind(Vmat, v)
  }
  
  B = tcrossprod(Rmat, Qmat) # PLS regression coefficients
  
  return(list(
    R = Rmat,
    coef = B
  ))
}

#' Principal centroid classification
#' x.old is a matrix with rows for subjects and columns corresponding to time points
#' y.old is a vector of 0 and 1
#' par is a matrix with 2 rows: 1st row for candidate number of components and 2nd one for penalty theta

PCC = function(y.old, x.old, y.new = NULL, x.new = NULL, timePts, pars){

  m = ncol(x.old)
  if (is.null(x.new)) expansion = expand.bspline(x.old, timePts)
  else expansion = expand.bspline(rbind(x.old, x.new), timePts)
  splineBasis = expansion$basis
  basismatrix = eval.basis(timePts, expansion$basis)
  W = inprod(splineBasis, splineBasis)
  D = inprod(splineBasis, splineBasis, Lfdobj1 = 2, Lfdobj2 = 2)
  
  Whalf = half(W)
  
  errPredVec = NULL
  
  index.0 = (y.old==0)
  index.1 = (y.old==1)
  n.0 = sum(index.0)
  n.1 = sum(index.1)
  n = length(y.old)
  
  mu = colMeans(x.old)
  mu.0 = colMeans(x.old[index.0,])
  mu.1 = colMeans(x.old[index.1,])
  
  C.old.diff = t(smooth.basis(
    y=t(sweep(x.old, 2, mu)), 
    argvals=timePts, 
    expansion$D2fdPar, 
    dfscale=1.2
  )$fd$coefs)
  C.old.diff.resp = t(smooth.basis(
    y=t(rbind(sweep(x.old[index.0,], 2, mu.0), sweep(x.old[index.1,], 2, mu.1))), 
    argvals=timePts, 
    expansion$D2fdPar, 
    dfscale=1.2
  )$fd$coefs)
  
  if (!is.null(x.new)) {
    C.new.diff.0 = t(smooth.basis(
      y=t(sweep(x.new, 2, mu.0)), 
      argvals=timePts, 
      expansion$D2fdPar, 
      dfscale=1.2
    )$fd$coefs)
    C.new.diff.1 = t(smooth.basis(
      y=t(sweep(x.new, 2, mu.1)), 
      argvals=timePts, 
      expansion$D2fdPar, 
      dfscale=1.2
    )$fd$coefs)
  }
  
  for (i in 1:ncol(pars)){
    p = pars[1, i]
    theta = pars[2, i]
    
    U = crossprod(C.old.diff.resp %*% W)/n
    G = W + theta*D
    G.halfinv = halfinv(G)
    V = G.halfinv %*% svd(crossprod(G.halfinv, U %*% G.halfinv), nu = p, nv = p)$u
    
    WV = W %*% V
    model = lm.fit(x = C.old.diff %*% WV, y = y.old - mean(y.old))
    
    if (!is.null(x.new)){
      D.pred = (C.new.diff.1 %*% WV %*% model$coef)^2/sum((Whalf %*% V %*% model$coef)^2) -
        (C.new.diff.0 %*% WV %*% model$coef)^2/sum((Whalf %*% V %*% model$coef)^2) +
        2*log(n.0/n.1)
      
      if (!is.null(y.new)) errPredVec = c(errPredVec, mean((y.new - (D.pred < 0))^2))
    } 
  }
  
  return(list(errPredVec = errPredVec))
}

#' Partial least centroid classification
#' x.old is a matrix with rows for subjects and columns corresponding to time points
#' y.old is a vector of 0 and 1
#' par is a vector of candidata number of components

PLCC = function(y.old, x.old, y.new = NULL, x.new = NULL, timePts, pars){
  
  m = ncol(x.old)
  
  if (is.null(x.new)) expansion = expand.bspline(x.old, timePts)
  else expansion = expand.bspline(rbind(x.old, x.new), timePts)
  splineBasis = expansion$basis
  basismatrix = eval.basis(timePts, expansion$basis)
  W = inprod(splineBasis, splineBasis)

  errPredVec = NULL
  Whalf = half(W)
  Whalfinv = halfinv(W)
  
  index.0 = (y.old==0)
  index.1 = (y.old==1)
  n.0 = sum(index.0)
  n.1 = sum(index.1)

  mu = colMeans(x.old)
  mu.0 = colMeans(x.old[index.0,])
  mu.1 = colMeans(x.old[index.1,])
  
  C.old.diff = t(smooth.basis(
    y=t(sweep(x.old, 2, mu)), 
    argvals=timePts, 
    expansion$D2fdPar, 
    dfscale=1.2
  )$fd$coefs)
  
  if (!is.null(x.new)) {
    C.new.diff.0 = t(smooth.basis(
      y=t(sweep(x.new, 2, mu.0)), 
      argvals=timePts, 
      expansion$D2fdPar, 
      dfscale=1.2
    )$fd$coefs)
    C.new.diff.1 = t(smooth.basis(
      y=t(sweep(x.new, 2, mu.1)), 
      argvals=timePts, 
      expansion$D2fdPar, 
      dfscale=1.2
    )$fd$coefs)
  }

  for (i in 1:length(pars)){
    p = pars[i]
    
    simpls.coef = simpls(y = matrix(y.old, ncol = 1), x = C.old.diff %*% Whalf, ncomp = p)$coef
    beta.coef = Whalfinv %*% simpls.coef
    beta.norm = as.numeric(t(beta.coef) %*% W %*% beta.coef)^.5
    
    if (!is.null(x.new)){
      D.pred = (C.new.diff.1 %*% half(W) %*% simpls.coef)^2 / beta.norm^2 -
        (C.new.diff.0 %*% half(W) %*% simpls.coef)^2 / beta.norm^2 +
        2*log(n.0/n.1)
      
      if (!is.null(y.new)) errPredVec[p] = mean(abs(y.new - (D.pred<0)) > 1e-1)
    } 
  }
  
  return(list(errPredVec = errPredVec))
}

#' CCC is implementing CCC with 1-dim functional x and 1-dim scalar y.
#' @param y is 1-dim 0-1 response.
#' @param x is 1-dim independent variable.
#' @param Alpha consists of elements in [0,1).
#' @param pMax is the initial max number of basis functions.

CCC = function(y.old, x.old, y.new = NULL, x.new = NULL, timePts, Alpha, pMax = 5, pMin = 1, random.search = T){

  m = ncol(x.old)
  n = nrow(x.old)
  y.old.cen = y.old - mean(y.old) # centering y.old
  
  index0 = (y.old==0)
  index1 = (y.old==1)
  n0 = sum(index0)
  n1 = sum(index1)
  
  mu = colMeans(x.old)
  mu0 = colMeans(x.old[index0,])
  mu1 = colMeans(x.old[index1,])
  
  if (is.null(x.new)) 
    expansion = expand.bspline(x.old, timePts)
  else 
    expansion = expand.bspline(rbind(x.old, x.new), timePts)

  splineBasis = expansion$basis
  basismatrix = eval.basis(timePts, expansion$basis)
  W = inprod(splineBasis, splineBasis)
  
  C.old.diff = t(smooth.basis(
    y=t(sweep(x.old, 2, mu)), 
    argvals=timePts, 
    expansion$D2fdPar, 
    dfscale=1.2
  )$fd$coefs)
  C.old.diff.0 = t(smooth.basis(
    y=t(sweep(x.old, 2, mu0)), 
    argvals = timePts, 
    expansion$D2fdPar, 
    dfscale=1.2
  )$fd$coefs)
  C.old.diff.1 = t(smooth.basis(
    y=t(sweep(x.old, 2, mu1)),
    argvals= timePts, 
    expansion$D2fdPar, 
    dfscale=1.2
  )$fd$coefs)

  if (!is.null(x.new)){
    C.new.diff.0 = t(smooth.basis(
      y=t(sweep(x.new, 2, mu0)), 
      argvals=timePts,
      expansion$D2fdPar, 
      dfscale=1.2
    )$fd$coefs)
    C.new.diff.1 = t(smooth.basis(
      y=t(sweep(x.new, 2, mu1)), 
      argvals=timePts, 
      expansion$D2fdPar, 
      dfscale=1.2
    )$fd$coefs)
  }

  # initiation
  betaCoef = array(NA, dim=c(length(Alpha), pMax, ncol(W))) 
  
  GCV.CCC = matrix(NA, ncol=pMax, nrow=length(Alpha))
  GCV.CCCL = matrix(NA, ncol=pMax, nrow=length(Alpha))
  GCV.CCCQ = matrix(NA, ncol=pMax, nrow=length(Alpha))
  GCV.LCC = matrix(NA, ncol=pMax, nrow=length(Alpha))
  GCV.QCC = matrix(NA, ncol=pMax, nrow=length(Alpha))
  
  errPred.CCC = matrix(NA, ncol=pMax, nrow=length(Alpha))
  errPred.CCCL = matrix(NA, ncol=pMax, nrow=length(Alpha))
  errPred.CCCQ = matrix(NA, ncol=pMax, nrow=length(Alpha))
  errPred.LCC = matrix(NA, ncol=pMax, nrow=length(Alpha))
  errPred.QCC = matrix(NA, ncol=pMax, nrow=length(Alpha))
  
  if (random.search) iMax = sample(pMin:pMax, size = length(Alpha), replace = T)
  else iMax = rep(pMax, times = length(Alpha))

  for(j in 1:length(Alpha)){
    
    alpha = Alpha[j]
    gam = alpha / (1-alpha)
    
    b = list() # initiate b
    G = list() # initiate G
    H = matrix(NA, nrow = n, ncol = iMax[j]) # initiate H
    wCoef = matrix(NA, nrow = ncol(W), ncol = iMax[j])
    
    Whalf = half(W)
    Whalfinv = halfinv(W)
    CWhalf = C.old.diff %*% Whalf
    PCWhalf = CWhalf
    r = qr(CWhalf)$rank
    V = list()

    for(i in 1:iMax[j]){
      
      if (i==1)
        P = diag(n)
      else
        P = P %*% (diag(n) - H[, i-1] %*% crossprod(H[, i-1])^{-1} %*% t(H[, i-1]))
      
      PCWhalf = P %*% PCWhalf
      
      if (i == 1){
        svdPCWhalf = svd(PCWhalf)
        V[[i]] = svdPCWhalf$v[, 1:r]
      }else
        V[[i]] = V[[1]]
      G[[i]] = PCWhalf %*% V[[i]]
      Gy = crossprod(G[[i]], y.old.cen)
      GG = crossprod(G[[i]])
      zeta = max(eigen(GG, only.values = T)$values)
  	  
      L = function(delta) {
        return(solve(GG + delta^-1 * zeta * diag(r)))
      }
  	  
      neglnQ = function(delta) {
        Ldelta = L(delta)
        result = -as.numeric(
            2 * log(abs(t(Gy) %*% Ldelta %*% Gy)) - 
            gam * log(abs(t(Gy) %*% Ldelta %*% Ldelta %*% Gy)) +
            (gam-1) * log(abs(t(Gy) %*% Ldelta %*% GG %*% Ldelta %*% Gy))
          )
        if (is.infinite(result)) return(1e10)
        else return(result)
      }
      
      deltaOpt = manymodes(neglnQ, ninitial = ninitial)
      LdeltaOpt = L(deltaOpt)
      b[[i]] = LdeltaOpt %*% Gy / as.numeric(crossprod(LdeltaOpt %*% Gy))^.5
      
      if (length(b)<i) break
      H[, i] = CWhalf %*% V[[i]] %*% b[[i]]
      wCoef[, i] = Whalfinv %*% V[[i]] %*% b[[i]]
      
      betaCoef[j ,i, ] = as.matrix(wCoef[, 1:i]) %*% lm.fit(x = as.matrix(H[, 1:i]), y = y.old.cen)$coefficients
      
      sigma0 = sd(C.old.diff.0[index0,] %*% W %*% betaCoef[j ,i,])
      sigma1 = sd(C.old.diff.1[index1,] %*% W %*% betaCoef[j ,i,])
      sigmaPool = (n-2)^-.5 * ((n0-1) * sigma0^2 + (n1-1) * sigma1^2)^.5
      Sigma0 = var(C.old.diff.0[index0,] %*% W %*% wCoef)
      Sigma1 = var(C.old.diff.1[index1,] %*% W %*% wCoef)
      SigmaPool = (n-2)^-1 * ((n0-1) * Sigma0 + (n1-1) * Sigma1)
      
      betaNorm = as.numeric(t(betaCoef[j ,i, ]) %*% W %*% betaCoef[j ,i, ])
      
      Dfit.CCC = betaNorm^-2 * (C.old.diff.1 %*% W %*% betaCoef[j ,i,])^2 - 
        betaNorm^-2 * (C.old.diff.0 %*% W %*% betaCoef[j ,i,])^2 + 
        2*log(n0/n1)	
      Dfit.CCCL = sigmaPool^-2 * (C.old.diff.1 %*% W %*% betaCoef[j ,i,])^2 - 
        sigmaPool^-2 * (C.old.diff.0 %*% W %*% betaCoef[j ,i,])^2 + 
        2*log(n0/n1)
      Dfit.CCCQ = sigma1^-2 * (C.old.diff.1 %*% W %*% betaCoef[j ,i,])^2 - 
        sigma0^-2 * (C.old.diff.0 %*% W %*% betaCoef[j ,i,])^2 + 
        2*log(n0 * sigma1 / n1 / sigma0)
      
      if (class(try(solve(SigmaPool),silent=T)) == "matrix"){
        Dfit.LCC = diag(C.old.diff.1 %*% W %*% wCoef %*% solve(SigmaPool) %*% t(C.old.diff.1 %*% W %*% wCoef)) -
          diag(C.old.diff.0 %*% W %*% wCoef %*% solve(SigmaPool) %*% t(C.old.diff.0 %*% W %*% wCoef)) +
          2*log(n0/n1)
      }else Dfit.LCC = sample(c(0, 1), size = length(y.old), replace = T, prob = c(n0/n, n1/n))
      
      if (class(try(solve(Sigma0), silent=T)) == "matrix" & class(try(solve(Sigma1), silent=T)) == "matrix"){
        Dfit.QCC = diag(C.old.diff.1 %*% W %*% wCoef %*% solve(Sigma1) %*% t(C.old.diff.1 %*% W %*% wCoef)) -
          diag(C.old.diff.0 %*% W %*% wCoef %*% solve(Sigma0) %*% t(C.old.diff.0 %*% W %*% wCoef)) +
          log(det(Sigma1)) - log(det(Sigma0)) + 2*log(n0/n1)
      }else Dfit.QCC = sample(c(0, 1), size = length(y.old), replace = T, prob = c(n0/n, n1/n))
              
      GCV.CCC[j, i] = sum(abs(y.old - (Dfit.CCC<0)) > 1e-1) / (n-i-2)^2
      GCV.CCCL[j, i] = sum(abs(y.old - (Dfit.CCCL<0)) > 1e-1) / (n-i-2)^2
      GCV.CCCQ[j, i] = sum(abs(y.old - (Dfit.CCCQ<0)) > 1e-1) / (n-i-2)^2
      GCV.LCC[j, i] = sum(abs(y.old - (Dfit.LCC<0)) > 1e-1) / (n-i-2)^2
      GCV.QCC[j, i] = sum(abs(y.old - (Dfit.QCC<0)) > 1e-1) / (n-i-2)^2
      
      if (!is.null(x.new)) {
        Dpred.CCC = betaNorm^-2 * (C.new.diff.1 %*% W %*% betaCoef[j ,i,])^2 - 
          betaNorm^-2 * (C.new.diff.0 %*% W %*% betaCoef[j ,i,])^2 + 
          2*log(n0/n1)	
        Dpred.CCCL = sigmaPool^-2 * (C.new.diff.1 %*% W %*% betaCoef[j ,i,])^2 - 
          sigmaPool^-2 * (C.new.diff.0 %*% W %*% betaCoef[j ,i,])^2 + 
          2*log(n0/n1)
        Dpred.CCCQ = sigma1^-2 * (C.new.diff.1 %*% W %*% betaCoef[j ,i,])^2 - 
          sigma0^-2 * (C.new.diff.0 %*% W %*% betaCoef[j ,i,])^2 + 
          2*log(n0 * sigma1 / n1 / sigma0)
        
        if (class(try(solve(SigmaPool), silent=T)) == "matrix"){
          Dpred.LCC = diag(C.new.diff.1 %*% W %*% wCoef %*% solve(SigmaPool) %*% t(C.new.diff.1 %*% W %*% wCoef)) -
            diag(C.new.diff.0 %*% W %*% wCoef %*% solve(SigmaPool) %*% t(C.new.diff.0 %*% W %*% wCoef)) +
            2*log(n0/n1)
        }else Dpred.LCC = sample(c(0, 1), size = length(y.new), replace = T, prob = c(n0/n, n1/n))
        
        if (class(try(solve(Sigma0), silent=T)) == "matrix" & class(try(solve(Sigma1), silent=T)) == "matrix"){
          Dpred.QCC = diag(C.new.diff.1 %*% W %*% wCoef %*% solve(Sigma1) %*% t(C.new.diff.1 %*% W %*% wCoef)) -
            diag(C.new.diff.0 %*% W %*% wCoef %*% solve(Sigma0) %*% t(C.new.diff.0 %*% W %*% wCoef)) +
            log(det(Sigma1)) - log(det(Sigma0)) + 2*log(n0/n1)
        }else Dpred.QCC = sample(c(0, 1), size = length(y.new), replace = T, prob = c(n0/n, n1/n))
        
    		if (!is.null(y.new)) {
    		  errPred.CCC[j, i] = mean((y.new - (Dpred.CCC<0))^2)
    		  errPred.CCCL[j, i] = mean((y.new - (Dpred.CCCL<0))^2)
    		  errPred.CCCQ[j, i] = mean((y.new - (Dpred.CCCQ<0))^2)
    		  errPred.LCC[j, i] = mean((y.new - (Dpred.LCC<0))^2)
    		  errPred.QCC[j, i] = mean((y.new - (Dpred.QCC<0))^2)
    		}
      }
    }
  }
  
  locOpt.CCC = which(GCV.CCC == min(GCV.CCC, na.rm=TRUE), arr.ind=TRUE)
  locOpt.CCCL = which(GCV.CCCL == min(GCV.CCCL, na.rm=TRUE), arr.ind=TRUE)
  locOpt.CCCQ = which(GCV.CCCQ == min(GCV.CCCQ, na.rm=TRUE), arr.ind=TRUE)
  locOpt.LCC = which(GCV.LCC == min(GCV.LCC, na.rm=TRUE), arr.ind=TRUE)
  locOpt.QCC = which(GCV.QCC == min(GCV.QCC, na.rm=TRUE), arr.ind=TRUE)
  
  parsOpt.CCC = NULL
  parsOpt.CCCL = NULL
  parsOpt.CCCQ = NULL
  parsOpt.LCC = NULL
  parsOpt.QCC = NULL

  betaCoefOpt.CCC = matrix(NA, ncol=nrow(locOpt.CCC), nrow=ncol(W))
  betaCoefOpt.CCCL = matrix(NA, ncol=nrow(locOpt.CCCL), nrow=ncol(W))
  betaCoefOpt.CCCQ = matrix(NA, ncol=nrow(locOpt.CCCQ), nrow=ncol(W))
  betaCoefOpt.LCC = matrix(NA, ncol=nrow(locOpt.LCC), nrow=ncol(W))
  betaCoefOpt.QCC = matrix(NA, ncol=nrow(locOpt.QCC), nrow=ncol(W))
  
  errPredOpt.CCC = NULL
  errPredOpt.CCCL = NULL
  errPredOpt.CCCQ = NULL
  errPredOpt.LCC = NULL
  errPredOpt.QCC = NULL

  for (i in 1:nrow(locOpt.CCC)){ # in case there are multi optimizers
    parsOpt.CCC = rbind(parsOpt.CCC, c(Alpha[locOpt.CCC[i,1]], locOpt.CCC[i,2]))
    betaCoefOpt.CCC[,i] = betaCoef[locOpt.CCC[i,1], locOpt.CCC[i,2], ]
  }
  colnames(parsOpt.CCC) = c("alpha", "ncomp")
  
  for (i in 1:nrow(locOpt.CCCL)){ # in case there are multi optimizers
    parsOpt.CCCL = rbind(parsOpt.CCCL, c(Alpha[locOpt.CCCL[i,1]], locOpt.CCCL[i,2]))
    betaCoefOpt.CCCL[,i] = betaCoef[locOpt.CCCL[i,1], locOpt.CCCL[i,2], ]
  }
  colnames(parsOpt.CCCL) = c("alpha", "ncomp")
  
  for (i in 1:nrow(locOpt.CCCQ)){ # in case there are multi optimizers
    parsOpt.CCCQ = rbind(parsOpt.CCCQ, c(Alpha[locOpt.CCCQ[i,1]], locOpt.CCCQ[i,2]))
    betaCoefOpt.CCCQ[,i] = betaCoef[locOpt.CCCQ[i,1], locOpt.CCCQ[i,2], ]
  }
  colnames(parsOpt.CCCQ) = c("alpha", "ncomp")
  
  for (i in 1:nrow(locOpt.LCC)){ # in case there are multi optimizers
    parsOpt.LCC = rbind(parsOpt.LCC, c(Alpha[locOpt.LCC[i,1]], locOpt.LCC[i,2]))
    betaCoefOpt.LCC[,i] = betaCoef[locOpt.LCC[i,1], locOpt.LCC[i,2], ]
  }
  colnames(parsOpt.LCC) = c("alpha", "ncomp")
  
  for (i in 1:nrow(locOpt.QCC)){ # in case there are multi optimizers
    parsOpt.QCC = rbind(parsOpt.QCC, c(Alpha[locOpt.QCC[i,1]], locOpt.QCC[i,2]))
    betaCoefOpt.QCC[,i] = betaCoef[locOpt.QCC[i,1], locOpt.QCC[i,2], ]
  }
  colnames(parsOpt.QCC) = c("alpha", "ncomp")
  
  if (!is.null(x.new) & !is.null(y.new)){
    for (i in 1:nrow(locOpt.CCC)){ # in case there are multi optimizers
      if (i==1) errPredOpt.CCC = errPred.CCC[locOpt.CCC[i,1], locOpt.CCC[i,2]]
	    else errPredOpt.CCC = c(errPredOpt.CCC, errPred.CCC[locOpt.CCC[i,1], locOpt.CCC[i,2]])
    }
	
    for (i in 1:nrow(locOpt.CCCL)){ # in case there are multi optimizers
      if (i==1) errPredOpt.CCCL = errPred.CCCL[locOpt.CCCL[i,1], locOpt.CCCL[i,2]]
      else errPredOpt.CCCL = c(errPredOpt.CCCL, errPred.CCCL[locOpt.CCCL[i,1], locOpt.CCCL[i,2]])
    }
    
    for (i in 1:nrow(locOpt.CCCQ)){ # in case there are multi optimizers
      if (i==1) errPredOpt.CCCQ = errPred.CCCQ[locOpt.CCCQ[i,1], locOpt.CCCQ[i,2]]
      else errPredOpt.CCCQ = c(errPredOpt.CCCQ, errPred.CCCQ[locOpt.CCCQ[i,1], locOpt.CCCQ[i,2]])
    }
    
    for (i in 1:nrow(locOpt.LCC)){ # in case there are multi optimizers
      if (i==1) errPredOpt.LCC = errPred.LCC[locOpt.LCC[i,1], locOpt.LCC[i,2]]
      else errPredOpt.LCC = c(errPredOpt.LCC, errPred.LCC[locOpt.LCC[i,1], locOpt.LCC[i,2]])
    }
    
    for (i in 1:nrow(locOpt.QCC)){ # in case there are multi optimizers
      if (i==1) errPredOpt.QCC = errPred.QCC[locOpt.QCC[i,1], locOpt.QCC[i,2]]
      else errPredOpt.QCC = c(errPredOpt.QCC, errPred.QCC[locOpt.QCC[i,1], locOpt.QCC[i,2]])
    }
  }
  
  return(list(GCV.CCC = GCV.CCC, GCV.CCCL = GCV.CCCL, GCV.CCCQ = GCV.CCCQ, GCV.LCC = GCV.LCC, GCV.QCC = GCV.QCC,
              D2fdPar = expansion$D2fdPar, 
              errPredOpt.CCC= errPredOpt.CCC, 
              errPredOpt.CCCL = errPredOpt.CCCL, errPredOpt.CCCQ= errPredOpt.CCCQ,
              errPredOpt.LCC= errPredOpt.LCC, errPredOpt.QCC= errPredOpt.QCC,
      			  parsOpt.CCC = parsOpt.CCC, 
      			  parsOpt.CCCL = parsOpt.CCCL, parsOpt.CCCQ = parsOpt.CCCQ,
      			  parsOpt.LCC = parsOpt.LCC, parsOpt.QCC = parsOpt.QCC,
      			  betaCoefOpt.CCC = betaCoefOpt.CCC, 
      			  betaCoefOpt.CCCL = betaCoefOpt.CCCL, betaCoefOpt.CCCQ = betaCoefOpt.CCCQ,
      			  betaCoefOpt.LCC = betaCoefOpt.LCC, betaCoefOpt.QCC = betaCoefOpt.QCC
      			  ))
}

#' integrate error rates

creatErrMat = function(errRateLst){
  nrowErrRateMat = min(lengths(errRateLst))
  errRateMat = NULL
  if (nrowErrRateMat > 0){
    for (i in 1:length(errRateLst)){
      errRateMat = cbind(errRateMat, 100 * errRateLst[[i]][1:nrowErrRateMat])
    }
    colnames(errRateMat) = names(errRateLst)
  }
  return(errRateMat)
}

#'  boxplot of error rate

boxplot.error.rate = function(nrowplot, error.rate){
  
  if (!("ggplot2" %in% rownames(installed.packages()))) install.packages("ggplot2")
  library(ggplot2)
  library(reshape2)
  
  melton = melt(
    data.frame(error.rate[1:nrowplot,], 
              Replication=1:nrowplot), 
              id.vars = "Replication")
  bplot = ggplot(melton, aes(x = variable, y = value, colour = variable)) + 
    geom_boxplot(outlier.shape=NA, width=.3) +
    coord_cartesian(ylim = c(0, 40)) +
    labs(x='', y='Misclassification %') + 
    ggtitle("") +
    theme_bw() +
    theme(legend.position = "none", 
          panel.border = element_blank(), 
          plot.title = element_text(hjust = .5),
          text = element_text(size = 25),
          axis.text.x = element_text(size = 25, angle = 60, hjust = .8)
          ) 
  
  plot(bplot)
  file = switch(simu+1,
                paste('bplot_',
                      switch(case.no, 
                             'TecatorProtein_',
                             'TecatorProteinDerive_',
                             'DTIcca5_'),  
                      nrowplot, 'repeats_', 
                      prop.train*100, 'proptrain_',
                      pMax, 'pMax.pdf', 
                      sep=''),
                paste('bplot_Simu_',
                      case.no, '_',
                      'pi0=', pi0*100, '_',
                      'rho=', rho, '_',
                      nrowplot, 'repeats_', 
                      prop.train*100, 'proptrain_',
                      pMax, 'pMax.pdf', 
                      sep='')
  )
  ggsave(file = file, 
         width = 6,
         height = 8,
         units = 'in',
         dpi = 300,
         path = "figure")
}

PlotMultiCurve = function(x, y){
  
  if (!("ggplot2" %in% rownames(installed.packages()))) install.packages("ggplot2")
  if (!("tidyr" %in% rownames(installed.packages()))) install.packages("tidyr")
  library(ggplot2)
  
  data.reform = data.frame(ID = 1:nrow(x), class = factor(y), time = x)
  data.reform = tidyr::gather(data.reform, time, value, -c(ID, class))
  data.reform$time = rep((1:ncol(x))/100, each = nrow(x))
  curvesplot.simu = 
    ggplot(data.reform, aes(x=time, y=value, colour=class, linetype=class)) +
    geom_line(aes(group=ID))
  curvesplot.simu +
    coord_cartesian(ylim = c(-250, 250)) +
    labs(x='t', y='X(t)') + 
    ggtitle("") +
    theme_bw()+
    theme(legend.position = "none", 
          panel.border = element_blank(), 
          plot.title = element_text(hjust = .5),
          text = element_text(size = 25)
          )
}