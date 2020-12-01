#' Put main.R & supportFun.R into the identical folder.
#' In the same folder, create two subfolders named figure' & 'Rimage', respectively.

if (!("rstudioapi" %in% rownames(installed.packages()))) install.packages("rstudioapi")
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("functions.R")

if (!("fda.usc" %in% rownames(installed.packages()))) install.packages("fda.usc")
if (!("e1071" %in% rownames(installed.packages()))) install.packages("e1071")
library(fda.usc)

#### global setting
set.seed(1)
RR = 200L # number of replication
prop.train = .8 # proportion of whole dataset left for training
ninitial = 1 # number of initial values used in the optimization for CCC
epsilon = 1e-9 # small non-zero values used in @function CCC 
nfold = 5 # numbe of folds for CV
simu = 1 # 1 for simulation, 0 for real data
case.no = 2
rho = 10 # signal-to-noise ratio (SNR) for simulation only

#### Simulation

### communal setting for simulation
TT = 100 # total number of time points
N = 200 # number of curves
pi0 = .8 # prior for population 0
J = 5 # number of eigenfunctions considered
distri = 3 # 1 for normal, 2 for (chi2(4)-4)/8^.5, 3 for exp(1)-1

## scene I: difference arises early with SNR rho

if (simu == 1 & case.no == 1) {

  timePts = (0:TT)/TT # time points evaluated for each curve
  
  Lambda.0 = c(200, 100, 1, .2, .1)
  Lambda.1 = Lambda.0
  J = length(Lambda.0)
  
  mu.0 = rep(0, J)
  mu.1 = c(1, rep(0, J-1)) * rho * Lambda.1^.5
  
  eigenfun.0 = orthoBasis(order = 1:J, denseGrid = timePts, type = 'shiftedLegendre', normalized = T)
  eigenfun.1 = eigenfun.0
  
  y = list()
  x = list()
  for (R in 1:RR){
    y[[R]] = matrix(rbinom(N, size = 1, 1-pi0))
    x[[R]] = matrix(NA, nrow=N, ncol=ncol(eigenfun.0))
    for (n in 1:N){
      Z = switch(distri, rnorm(J), (rgamma(J, shape = 2, rate = .5) - 4)/sqrt(8), rexp(J)-1)
      if (y[[R]][n]) x[[R]][n, ] = matrix(sqrt(Lambda.1) * Z + mu.1, nrow = 1) %*% eigenfun.1
      else x[[R]][n, ] = matrix(sqrt(Lambda.0) * Z + mu.0, nrow = 1) %*% eigenfun.0
    }
  }
  PlotMultiCurve(x[[1]], y[[1]])
}

## scene II: difference arises late with SNR rho

if (simu == 1 & case.no == 2){
  timePts = (0:TT)/TT # time points evaluated for each curve
  Lambda.0 = c(200, 100, 1, .2, .1)
  Lambda.1 = rev(Lambda.0)
  names(Lambda.0) = names(Lambda.1) = c('1', '2', '3', '4', '5')
  J = length(Lambda.0)
  mu.0 = rep(0, J)
  mu.1 = c(0, 0, 1, 0, 0) * rho * Lambda.1^.5
  
  eigenfun.0 = orthoBasis(order = 1:J, denseGrid = timePts, type = 'shiftedLegendre', normalized = T)
  eigenfun.1 = eigenfun.0
  
  y = list()
  x = list()
  for (R in 1:RR){
    y[[R]] = matrix(rbinom(N, 1, 1-pi0))
    x[[R]] = matrix(NA, nrow = N, ncol = ncol(eigenfun.0))
    for (n in 1:N){
      Z = switch(distri, rnorm(J), (rgamma(J, shape = 2, rate = .5)-4)/sqrt(8), rexp(J)-1)
      if (y[[R]][n]) x[[R]][n, ] = matrix(sqrt(Lambda.1)*Z + mu.1, nrow=1) %*% eigenfun.1
      else x[[R]][n, ] = matrix(sqrt(Lambda.0)*Z + mu.0, nrow=1) %*% eigenfun.0
    }
  }

  PlotMultiCurve(x[[1]], y[[1]])
}

#### Real data

## Tecator

if (simu == 0 & case.no == 1){
  meat = t(matrix(t(as.matrix(read.table("data/meat.txt"))), ncol=240)); x = meat[, 1:100]
  y = matrix(as.numeric(meat[, 125]<16)) #protein
  timePts = seq(from = 0, to = 1, length.out = ncol(x))
}

if (simu == 0 & case.no == 2){
  meat = t(matrix(t(as.matrix(read.table("data/meat.txt"))), ncol=240)); x = meat[, 1:100]
  y = matrix(as.numeric(meat[, 125]<16)) #protein
  timePts = seq(from = 0, to = 1, length.out = ncol(x))
  expansion = expand.bspline(x, timePts = timePts)
  x = t(fda::eval.fd(evalarg = timePts, fdobj = expansion$fdObj, Lfdobj = 2))
}

## DTI data from R package refund

if (simu == 0 & case.no == 3){
  if (!("refund" %in% rownames(installed.packages()))) install.packages("refund")
  attach(refund::DTI)
  y = matrix(as.numeric(case))
  x = impute.curve(cca)
  detach(refund::DTI)
  timePts = seq(from = 0, to = 1, length.out = ncol(x))
}

#######################################################################
#######################################################################
#######################################################################

# Dangerous! clear existing result
errRate.CCCL.rs = NULL; errRate.CCCQ.rs = NULL; alpha.CCCQ.rs = NULL; run.time.rs = NULL
errRate.CCCL.fs = NULL; errRate.CCCQ.fs = NULL; alpha.CCCQ.fs = NULL; run.time.fs = NULL
errRate.PLCC = NULL; run.time.PLCC = NULL
errRate.PCC = NULL; run.time.PCC = NULL
errRate.logit = NULL; run.time.logit = NULL
errRate.nb = NULL; run.time.nb = NULL

# Check current progress
check.CCCL.rs = ifelse(is.null(errRate.CCCL.rs), 0, length(errRate.CCCL.rs))
check.CCCQ.rs = ifelse(is.null(errRate.CCCQ.rs), 0, length(errRate.CCCQ.rs))
check.CCCL.fs = ifelse(is.null(errRate.CCCL.fs), 0, length(errRate.CCCL.fs))
check.CCCQ.fs = ifelse(is.null(errRate.CCCQ.fs), 0, length(errRate.CCCQ.fs))
check.PLCC = ifelse(is.null(errRate.PLCC), 0, length(errRate.PLCC))
check.PCC = ifelse(is.null(errRate.PCC), 0, length(errRate.PCC))
check.logit = ifelse(is.null(errRate.logit), 0, length(errRate.logit))
check.nb = ifelse(is.null(errRate.nb), 0, length(errRate.nb))

################ Replica ##################

if (!simu) pMax = pUpper.compu(x, y, timePts, .99) # max number of components for CCC, PCC and PLCC

for (R in 1:RR){
  
  # candidate pools for parameters to be tuned

  if (simu) pMax = J
  pMin = 1 # min number of components
  Alpha = c(seq(from = 0, to = .9, length.out = 10), .999) # supervision parameters for CCC
  Theta = c(0) # penalty for smoothness for PCC
  parsPCC = t(expand.grid(pMin:pMax, Theta)) # combinations of candidate parameters for PCC
  parsPLCC = pMin:pMax # candidate parameters for PLCC
  pars = list(parsPCC, parsPLCC)
  names(pars) = c('PCC', 'PLCC')

  if (!simu) {
    samp.idx = sample(1:nrow(x), round(nrow(x) * prop.train))
    x.old = x[samp.idx, ]
    y.old = matrix(y[samp.idx,])
    x.new = x[-samp.idx, ]
    y.new = matrix(y[-samp.idx,])
  }else{
    samp.idx = sample(1:nrow(x[[R]]), round(nrow(x[[R]]) * prop.train))
    x.old = x[[R]][samp.idx, ]
    y.old = matrix(y[[R]][samp.idx,])
    x.new = x[[R]][-samp.idx, ]
    y.new = matrix(y[[R]][-samp.idx,])
  }
  
  if (R > min(check.CCCL.rs, check.CCCQ.rs, check.CCCL.fs, check.CCCQ.fs)){
    ptm0 = proc.time()[3]
    result.rs = CCC(y.old, x.old, y.new, x.new, timePts, Alpha, pMax, pMin, random.search = T)
    errRate.CCCL.rs = c(errRate.CCCL.rs, mean(result.rs$errPredOpt.CCCL, na.rm = T))
    errRate.CCCQ.rs = c(errRate.CCCQ.rs, mean(result.rs$errPredOpt.CCCQ, na.rm = T))
    alpha.CCCQ.rs = c(alpha.CCCQ.rs, result.rs$parsOpt.CCCQ[, 1])
    run.time.rs[R] = proc.time()[3] - ptm0
    
    ptm0 = proc.time()[3]
    result.fs = CCC(y.old, x.old, y.new, x.new, timePts, Alpha, pMax, pMin, random.search = F)
    errRate.CCCL.fs = c(errRate.CCCL.fs, mean(result.fs$errPredOpt.CCCL, na.rm = T))
    errRate.CCCQ.fs = c(errRate.CCCQ.fs, mean(result.fs$errPredOpt.CCCQ, na.rm = T))
    alpha.CCCQ.fs = c(alpha.CCCQ.fs, result.fs$parsOpt.CCCQ[, 1])
    run.time.fs[R] = proc.time()[3] - ptm0
  }
  
  if (R > check.PLCC) {
    ptm0 = proc.time()[3]
    result = crossVali(y.old, x.old, y.new, x.new, timePts, nfold = 5, method = c('PLCC'), pars)
    errRate.PLCC = c(errRate.PLCC, mean(result$errPred$PLCC, na.rm = T))
    run.time.PLCC[R] = proc.time()[3] - ptm0
  }
  
  if (R > check.PCC) {
    ptm0 = proc.time()[3]
    result = crossVali(y.old, x.old, y.new, x.new, timePts, nfold = 5, method = c('PCC'), pars)
    errRate.PCC = c(errRate.PCC, mean(result$errPred$PCC, na.rm = T))
    run.time.PCC[R] = proc.time()[3] - ptm0
  }
  
  if (R > min(check.logit)) {
    ptm0 = proc.time()[3]
    x.old.fdobj = fdata(x.old, argvals = timePts)
    x.new.fdobj = fdata(x.new, argvals = timePts)
    out.logit = classif.glm(
      y ~ x, 
      data = list("df" = data.frame(y = factor(y.old)), "x" = x.old.fdobj)
    )
    pred.logit = predict(out.logit, list('x' = x.new.fdobj))
    errRate.logit = c(errRate.logit, mean(pred.logit != y.new))
    run.time.logit[R] = proc.time()[3] - ptm0
  }
  
  if (R > min(check.nb)) {
    ptm0 = proc.time()[3]
    x.old.fdobj = fdata(x.old, argvals = timePts)
    x.new.fdobj = fdata(x.new, argvals = timePts)
    out.nb = classif.naiveBayes(
      y ~ x, 
      data = list("df" = data.frame(y = factor(y.old)), "x" = x.old.fdobj)
    )
    pred.nb = predict(out.nb, list("x"=x.new.fdobj))
    errRate.nb = c(errRate.nb, mean(pred.nb != y.new))
    run.time.nb[R] = proc.time()[3] - ptm0
  }
  
  file = switch(simu+1,
                paste('Rimage/',
                      switch(case.no, 
                             'TecatorProtein_',
                             'TecatorProteinDerive_',
                             'DTIcca5_'),  
                      RR, 'repeats_', 
                      pMax, 'pMax.RData', 
                      sep=''),
                paste('Rimage/Simu_',
                      case.no, '_',
                      'pi0=', pi0*100, '_',
                      'rho=', rho, '_',
                      RR, 'repeats_',
                      pMax, 'pMax.RData', 
                      sep='')
  )
  
  if (R %% 40 == 0) save.image(file = file)
  if (R %% 40 == 0) cat(R,'\n')
  else cat(R)
}

#' boxplots

errRateLst = list(
  errRate.CCCL.rs
  ,errRate.CCCQ.rs
  # ,errRate.CCCL.fs
  # ,errRate.CCCQ.fs
  ,errRate.PLCC
  ,errRate.PCC
  ,errRate.logit
  ,errRate.nb
)
names(errRateLst) = c(
  'CCC-L'
  ,'CCC-Q'
  # ,'CCC-L.FG'
  # ,'CCC-Q.FG'
  ,'PLCC'
  ,'PCC'
  ,'Logit'
  ,'Naive.Bayes'
)

errMat = creatErrMat(errRateLst)
round(apply(errMat, 2, sd, na.rm = T), digits = 2)
round(apply(errMat, 2, mean, na.rm = T), digits = 2)
boxplot.error.rate(RR, errMat)
