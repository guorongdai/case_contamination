.libPaths("/home/guorongdai/RLibrary")

suppressMessages({
  
  if (!require(np)) install.packages("/home/grad/rondai/ssl/np_0.60-10.tar.gz", 
                                     repos = NULL, type = .Platform$pkgType)
  
  if (!require(glmnet)) install.packages('glmnet', repos =
                                           'https://cran.revolutionanalytics.com/')
  
  
  library(MASS)
  
  library(np)
  
  library(glmnet)
})


#################################################################################
# A solver for the estimating equation in Dai et al..
# x is the covariate matrix of the association model.
# w is the covariate matrix of the phenotyping model.
# d = 0, 1 of 2 means an observation is a false case, real case or control.
# s = 0 or 1 means an observation is in the control pool or the case pool.
# r = 1 means the true status of a case is known; r = 0 means the true status of a case is unknown.
# k is the number of folds in the cross validation for self fitting.
# k = NULL means leave-one-out cross validation.
# bset is the set of candidate bandwidth.
# nd is the number of directions in dimension reduction.
# dr = T means dimension reduction should be conducted.
# rseed is the random seed for the cross validation.
# sp = T means the debiased semi-parametric estimator should be implemented.
# sp = F means the debaised parametric estimator should be implemented.
# intercept1 = T means an intercept term should be added to x. 
# intercept2 = T means an intercept term should be added to w. 
# fi = T means full imputation; fi = F means partial imputation.
# eps is tolerence of the estimation. maxout is the maximum number of iterations.
# intval is the initial values. The default is a vector of zero.
# The output is a list including at most five items:
# 1. $ theta0 (if included) : the intercept of the association model;
# 2. $ theta: the slope of the association model;
# 3. $ alpha0 (if included): the intercept of the phenotyping model;
# 4. $ alpha: the slope of the phenotyping model.
# 5. $ covariance: the asymptotic variance-covariance matrix of 
#                  (theta0, theta ^ T, alpha0, alpha ^ T) ^ T.

solver = function(x, w, d, s, r, k = NULL, nd = NULL, dr = F, drmethod = NULL, sp = F, intercept1 = T,
                  intercept2 = T, fi = T, eps = 1e-8, maxout = 20, intval = NULL)
{
  
  dyn.load("ccs.so")
  
  if (is.vector(x)) x = matrix(x, ncol = 1)
  
  if (is.vector(w)) w = matrix(w, ncol = 1)
  
  n = as.integer(nrow(x))
  
  
  w1 = w[(s == 1) & (r == 1), ]
  w2 = w[(s == 1) & (r == 0), ]
  d1 = d[(s == 1) & (r == 1)]
  
  u = numeric(n)
  
  v = numeric(n)
  
  if(sp == T)
  {
    
    trans = drtrans(w1, d1, w2, drmethod, nd, dr)
    w1.trans = trans $ x1
    
    np.bw = npcdensbw(xdat = w1.trans, ydat = as.factor(d1), nmulti = 1)
    bw.opt = np.bw $ xbw # optimal bandwidth for full data
    
    
    n1 = nrow(w1)
    n2 = nrow(w2)
    nk = floor(n1 / k)
    
    
    
    if(is.vector(w1.trans))
    {
      
      p = 1
      
    } else {
      
      p = ncol(w1.trans) # dimensionality (after dimension reduction if applied)
      
    }
    
    
    h = bw.opt * ( ( n1 / (n1 - nk) ) ^ ( 1 / (4 + p) ) ) 
    
    offset1 = numeric(n1)
    offset2 = numeric(n2)
    
    for(i in 1 : (k - 1))
    {
      
      ind = (nk * (i - 1) + 1) : (nk * i)
      ind1 = 1 : length(ind)
      
      ww1 = w1[-ind, ]
      dd1 = d1[-ind]
      
      ww2 = rbind(w1[ind, ], w2)
      
      trans.ds = drtrans(ww1, dd1, ww2, drmethod, nd, dr)
      ww1.trans = trans.ds $ x1
      ww2.trans = trans.ds $ x2
      
      np.reg = npreg(txdat = ww1.trans, tydat = dd1, bws = h, exdat = ww2.trans)
      fit = np.reg $ mean
      
      offset1[ind] = fit[ind1]
      offset2 = offset2 + fit[-ind1] / k
      
    }
    
    ind = (nk * (k - 1) + 1) : n1
    ind1 = 1 : length(ind)
    
    ww1 = w1[-ind, ]
    dd1 = d1[-ind]
    
    ww2 = rbind(w1[ind, ], w2)
    
    trans.ds = drtrans(ww1, dd1, ww2, drmethod, nd, dr)
    ww1.trans = trans.ds $ x1
    ww2.trans = trans.ds $ x2
    
    np.reg = npreg(txdat = ww1.trans, tydat = dd1, bws = h, exdat = ww2.trans)
    fit = np.reg $ mean
    
    offset1[ind] = fit[ind1]
    offset2 = offset2 + fit[-ind1] / k
    
    
    u[(s == 1) & (r == 0)] = hinv( offset2 )
    u[(s == 1) & (r == 1)] = hinv( offset1 )
    u[u == Inf] = 10
    u[u == -Inf] = -10
    
    v[(s == 1) & (r == 1)] = hinv( offset1 )
    v[v == Inf] = 10
    v[v == -Inf] = -10
    
  }
  
  
  
  if (intercept1) x = cbind(rep(1, n), x)
  
  if (intercept2) w = cbind(rep(1, n), w)
  
  r = as.double(r)
  s = as.double(s)
  d = as.double(d)
  u = as.double(u)
  v = as.double(v)
  
  p1 = as.integer(ncol(x))
  p2 = as.integer(ncol(w))
  
  x = matrix(as.double(x), nrow = n, ncol = p1)
  w = matrix(as.double(w), nrow = n, ncol = p2)
  
  p = as.integer(p1 + p2)
  
  par = as.double(numeric(p))
  cov = matrix(as.double(0), p, p)
  if( is.null(intval) ) intval = numeric(p)
  intval = as.double(intval)
  
  fi = as.integer(fi)
  
  eps = as.double(eps)
  maxout = as.integer(maxout)
  
  
  out = .Fortran("eesolver", par, cov, intval, r, s, d, u, 
                 v, x, w, n, p1, p2, p, fi, eps, maxout)
  
  solution = out[[1]]
  abetainv = out[[2]]
  
  ee.m = matrix(0, n, p)
  
  ######################## 
  # asymptotic covariance matrix
  theta = solution[1 : p1]
  alpha = solution[-(1 : p1)]
  
  loc = which(s * r == 1) 
  nv = length(loc)
  
  ee = matrix(0, nv, p)
  
  xx = x[loc, ]
  ww = w[loc, ]
  dd = d[loc]
  
  for (i in 1 : nv)
  {
    
    ee[i, 1 : p1] = hf(sum(ww[i, ] * alpha)) * bh(sum(xx[i, ] * theta)) * xx[i, ]
    ee[i, -(1 : p1)] = ( dd[i] - hf(sum(ww[i, ] * alpha)) ) * bh(sum(xx[i, ] * theta)) * ww[i, ]
    
  }
  
  ee.m[1 : nv, ] = ee
  
  covl = cov(ee)
  
  
  
  loc = which(s * (1 - r) == 1) 
  nn1 = length(loc)
  
  ee = matrix(0, nn1, p)
  
  xx = x[loc, ]
  ww = w[loc, ]
  dd = d[loc]
  
  for (i in 1 : nn1) ee[i, 1 : p1] = hf(sum(ww[i, ] * alpha)) * bh(sum(xx[i, ] * theta)) * xx[i, ]
  
  ee.m[(nv + 1) : (nv + nn1), ] = ee
  
  cov1 = cov(ee)
  
  
  loc = which(s == 0) 
  nn0 = length(loc)
  
  ee = matrix(0, nn0, p)
  
  xx = x[loc, ]
  ww = w[loc, ]
  dd = d[loc]
  rr = r[loc]
  
  for (i in 1 : nn0) ee[i, 1 : p1] = hf(sum(xx[i, ] * theta)) * xx[i, ]
  
  ee.m[(nv + nn1 + 1) : n, ] = ee
  
  cov0 = cov(ee)
  
  
  bbeta = nv * covl + nn1 * cov1 + nn0 * cov0
  
  varcov = abetainv %*% bbeta %*% t(abetainv) 
  
  ########################
  
  if (intercept1)
  {
    
    if (intercept2) 
    {
      
      res = list("theta0" = solution[1], "theta" = solution[2 : p1], 
                 "alpha0" = solution[p1 + 1], "alpha" = solution[(p1 + 2) : p])
      
    } else {
      
      res = list("theta0" = solution[1], "theta" = solution[2 : p1], 
                 "alpha" = solution[(p1 + 1) : p])
      
    }
    
  } else {
    
    if (intercept2) {
      
      res = list("theta" = solution[1 : p1], "alpha0" = solution[p1 + 1], 
                 "alpha" = solution[(p1 + 2) : p])
      
    } else {
      
      res = list("theta" = solution[1 : p1], "alpha" = solution[(p1+1) : p])
      
    }
    
  }
  
  res = c(res, list("covariance" = varcov, "abetainv" = abetainv, "ee.m" = ee.m))
  
  return(res)
  
}
#################################################################################

#################################################################################
# calculate the covariance matrix of PIEE
piee.cov = function(theta, alpha, x, w, r, s, d)
{
  
  x = addone(x)
  w = addone(w)
  p1 = ncol(x)
  p2 = ncol(w)
  p = p1 + p2
  n = nrow(x)
  
  ee.m = matrix(0, n, p)
  
  ######################## 
  # asymptotic covariance matrix
  loc = which(s * r == 1) 
  nv = length(loc)
  
  ee = matrix(0, nv, p)
  
  xx = x[loc, ]
  ww = w[loc, ]
  dd = d[loc]
  
  for (i in 1 : nv)
  {
    
    ee[i, 1 : p1] = hf(sum(ww[i, ] * alpha)) * bh(sum(xx[i, ] * theta)) * xx[i, ]
    ee[i, -(1 : p1)] = ( dd[i] - hf(sum(ww[i, ] * alpha)) ) * ww[i, ]
    
  }
  
  ee.m[1 : nv, ] = ee
  
  covl = cov(ee)
  
  
  
  loc = which(s * (1 - r) == 1) 
  nn1 = length(loc)
  
  ee = matrix(0, nn1, p)
  
  xx = x[loc, ]
  ww = w[loc, ]
  dd = d[loc]
  
  for (i in 1 : nn1) ee[i, 1 : p1] = hf(sum(ww[i, ] * alpha)) * bh(sum(xx[i, ] * theta)) * xx[i, ]
  
  ee.m[(nv + 1) : (nv + nn1), ] = ee
  
  cov1 = cov(ee)
  
  
  loc = which(s == 0) 
  nn0 = length(loc)
  
  ee = matrix(0, nn0, p)
  
  xx = x[loc, ]
  ww = w[loc, ]
  dd = d[loc]
  rr = r[loc]
  
  for (i in 1 : nn0) ee[i, 1 : p1] = hf(sum(xx[i, ] * theta)) * xx[i, ]
  
  ee.m[(nv + nn1 + 1) : n, ] = ee
  
  cov0 = cov(ee)
  
  
  bbeta = nv * covl + nn1 * cov1 + nn0 * cov0
  
  omega = matrix(0, p, p)
  
  delta = nv / n
  
  for(i in 1 : n)
  {
    si = s[i]
    ri = r[i]
    xi = x[i, ]
    di = d[i]
    
    omega[1 : p1, 1 : p1] = omega[1 : p1, 1 : p1] + ( (si * hf(sum(alpha * xi)) * dbh(sum(theta * xi)) - (1 - si) * dh(sum(theta * xi))) ) * (xi %*% t(xi))
    
    omega[1 : p1, -(1 : p1)] = omega[1 : p1, -(1 : p1)] + si * dh(sum(alpha * xi)) * bh(sum(theta * xi)) * (xi %*% t(xi)) 
    
    omega[-(1 : p1), -(1 : p1)] = omega[-(1 : p1), -(1 : p1)] - ri * si * dh(sum(alpha * xi)) * (xi %*% t(xi))
    
  }
  
  abeta = omega
  abetainv = solve(abeta)
  
  varcov = abetainv %*% bbeta %*% t(abetainv) 
  
  return(varcov)
  
}
#################################################################################

#################################################################################
# supervised dimension reduction for semi-supervised data
# x is the covariate matrix of the labeled set.
# y is the response vector of the labeled set.
# xnew is the covariate matrix of the unlabeled set.
# drmethod is the method for dimension reduction.
# drmethod = "linear" means using linear regression.
# drmethod = "logistic" means using logstic regression.
# drmethod = "save" means using SAVE on I(Y < theta).
# drmethod = "sir" means using SIR.
# drmethod = "save.con" means using SAVE on Y.
# drmethod = "plinear" means using penalized linear regression.
# drmethod = "plogistic" means using penalized logistic regression.
# drmethod = "psir" means using penalized SIR.
# r is the number of directions in dimension reduction.
# dr = T means applying dimension reduction via the semi-supervised SIR.
# There are two outputs: $x1 is the transformation of the labeled covariate matrix;
# $x2 is the transformation of the unlabaled covariate matrix.
drtrans = function(x, y, xnew, drmethod = NULL, r = NULL, dr = F)
{
  
  if(dr == T)
  {
    
    
    if(drmethod == "logistic")
    {
      
      slm = glmnet(x, y, family = "binomial",lambda = 0) 
      
      
      x1 = predict(slm, newx = x, type = "link")
      x2 = predict(slm, newx = xnew, type = "link")
      
    }
    
    if(drmethod == "save")
    {
      
      tm = drsave(x, y, r)
      
      x1 = x %*% tm
      x2 = xnew %*% tm
      
    }
    
    
    
    if(drmethod == "plogistic")
    {
      
      slm = cv.glmnet(x, y, family = "binomial") # 10-fold cv to choose lambda
      
      x1 = predict(slm, newx = x, s = "lambda.min", type = "link")
      x2 = predict(slm, newx = xnew, s = "lambda.min", type = "link")
      
    }
    
    
  } else {
    
    x1 = x
    x2 = xnew
    
  }
  
  res = list("x1" = x1, "x2" = x2)
  
  return(res)
  
}
#################################################################################





#################################################################################
# dimension reduction on a binary response via SAVE
# x is the covariate matrix. 
# y is the response vector. 
# r is the number of directions.
# The output is a p * r transformation matrix whose columns are the directions 
# where p is column number of x.
drsave = function(x, y, r)
{
  
  suppressMessages({
    if (!require(expm)) install.packages('expm', repos=
                                           'https://cran.revolutionanalytics.com/')
    
    library(expm)
  })
  
  n = nrow(x)
  p = ncol(x)
  
  mu = apply(x, 2, mean)
  sigma = t(x) %*% x / (n-1) - mu %*% t(mu) * n / (n-1)
  sqrtinv = solve(sqrtm(sigma))
  
  z = (x - matrix(rep(mu, n), nrow = n, byrow = T)) %*% sqrtinv
  
  f = mean(y)
  
  z0 = z[y == 0, ]
  z1 = z[y == 1, ]
  n0 = nrow(z0)
  n1 = nrow(z1)
  
  mu0 = apply(z0, 2, mean)
  sigma0 = t(z0) %*% z0 / (n0-1) - mu0 %*% t(mu0) * n0 / (n0-1)
  mm0 = diag(p) - sigma0
  mu1 = apply(z1, 2, mean)
  sigma1 = t(z1) %*% z1 / (n1 - 1) - mu1 %*% t(mu1) * n1 / (n1-1)
  mm1 = diag(p) - sigma1
  
  omega = (1 - f) * (mm0 %*% mm0) + f * (mm1 %*% mm1)
  
  decom = eigen(omega)
  location = order(abs(decom $ values), decreasing = T)[1 : r]
  
  m = as.matrix(decom $ vectors[, location])
  
  res = sqrtinv %*% m
  
  return(res)
}
#################################################################################


#################################################################################
# generate a random sample from a logistic model
# n is the sample size. p is the dimensionality of the covariate.
# beta is the slop. beta0 is the intercept.
# dist is the type of distribution. 
# dist = 1 means a centered uniform distribution. 
# dist = 2 menas a centered normal distribution.
# bd is the upper bound of the uniform distribuiton.
# sigma is the standard deviation of the normal distribution.
# rseed is the random seed.
# The output is a list, whose first item is the covariatex 
# and the second item is the response y.

generator = function(n, p, beta, beta0, dist, bd = NULL, sigma = NULL, rseed)
{
  
  set.seed(rseed)
  
  if (dist == 1)
  {
    
    x = matrix(runif(n * p, -bd, bd), nrow = n, ncol = p)
    
  } else if (dist ==2 ) {
    
    x.raw = rnorm(2 * n * p, 0, sigma)
    x.raw = x.raw[abs(x.raw) < 5]
    x = matrix(x.raw[1 : (n * p)], nrow = n, ncol = p)
    
  } else {
    
    x = matrix(rnorm(n * p, 0, sigma), nrow = n, ncol = p)
    
  }
  
  prob = hf( (x %*% beta)[, 1] + beta0 )
  
  y = rbinom(n, 1, prob) 
  
  
  res = list("x" = x, "y" = y)
  
  return(res)
  
}

#################################################################################


#################################################################################
# logistic function
hf = function(x) 1 / (1 + exp(-x))

# 1 - logistic function
bh = function(x) 1 / (1 + exp(x))

# derivarive of hf
dh = function(x) exp(x) / ((1 + exp(x)) ^ 2)

# derivarive of bh
dbh = function(x) - exp(x) / ((1 + exp(x)) ^ 2)

# inverse logistic function
hinv = function(x) log(x / (1 - x))


# add a column of 1 to the left of a matrix
addone = function(x)
{
  
  n = nrow(x)
  res = cbind(rep(1, n), x)
  
  return(res)
  
}
#################################################################################

#################################################################################
# add a column of ones to the left of the input matrix
addone = function(x) cbind(rep(1, nrow(x)), x)
#################################################################################


adap = function(x, w, d, s, r, k = NULL, bset = NULL, nd = NULL, dr = F, rseed = NULL, 
                sp = F, intercept1 = T, intercept2 = T, fi = T, eps = 1e-8, maxout = 20, 
                intval = NULL)
{
  
  n = length(s)
  nv = sum(s * r)
  
  tau = sum(s) / n
  delta = nv / n
  
  validset = which(r == 1 & d != 0) # controls and real cases in the validation set
  wei = rep(1, n)
  wei[r == 1 & d == 1] = tau / delta
  
  if(intercept1 == 1)
  {
    tvt = glm(s ~ x, family = binomial, weights = wei, subset = validset) $ coefficients
  } else{
    tvt = glm(s ~ -1 + x, family = binomial, weights = wei, subset = validset) $ coefficients
  }
  
  
  sp = solver(x = x, w = w, d = d, s = s, r = r, sp = sp, intercept1 = intercept1,
              intercept2 = intercept2, fi = fi, eps = eps, maxout = maxout, intval = intval)
  
  if(intercept1 == 1)
  {
    hvt = c(sp $ theta0, sp $ theta)
  } else{
    hvt = sp $ theta
  }
  
  abetainv = sp $ abetainv
  ee.m = sp $ ee.m
  
  theta = tvt ####
  
  
  
  if (is.vector(x)) x = matrix(x, ncol = 1)
  if (is.vector(w)) w = matrix(w, ncol = 1)
  
  if (intercept1) x = addone(x)
  if (intercept2) w = addone(x)
  
  p1 = ncol(x)
  p2 = ncol(w)
  
  var.uiee = sp $ covariance[1 : p1, 1 : p1] # variance of UIEE
  
  
  gg = diag(tau * r * s * d * dbh((x %*% theta)[, 1]) - delta * (1 - s) * dh((x %*% theta)[, 1]))
  
  g = t(x) %*% gg %*% x 
  ginv = solve(g)
  
  
  ainv = abetainv[1 : p1, ]
  
  
  
  ########################################
  # variance of label-only and covariance between UIEE and label-only
  loc = which(s * r == 1) 
  
  ee = matrix(0, nv, p1)
  
  xx = x[loc, ]
  dd = d[loc]
  
  for (i in 1 : nv) ee[i, ] = tau * dd[i] * xx[i, ] * bh(sum(xx[i, ] * theta))
  
  ee.label = ee %*% t(ginv)
  varl = cov(ee.label)
  
  ee.uiee = ee.m[1 : nv, ] %*% t(ainv)
  covl = cov(ee.label, ee.uiee)
  
  
  loc = which(s == 0) 
  n0 = length(loc)
  
  ee = matrix(0, n0, p1)
  
  xx = x[loc, ]
  
  for (i in 1 : n0) ee[i, 1 : p1] = delta * xx[i, ] * hf(sum(xx[i, ] * theta))
  
  ee.label = ee %*% t(ginv)
  var0 = cov(ee.label)
  
  ee.uiee = ee.m[(n - n0 + 1) : n, ] %*% t(ainv)
  cov0 = cov(ee.label, ee.uiee)
  
  var.label = varl * nv + var0 * n0 # varaince of label-only
  cov.label.uiee = covl * nv + cov0 * n0 # covariance between label-only and UIEE
  
  ########################################
  
  var.label = diag(var.label)
  var.uiee = diag(var.uiee)
  cov.label.uiee = diag(cov.label.uiee)
  
  adapw = (var.label - cov.label.uiee) / (var.uiee + var.label - 2 * cov.label.uiee)
  
  
  adaptive = tvt + adapw * (hvt - tvt)
  var.adaptive = var.label + adapw ^ 2 * (var.uiee + var.label - 2 * cov.label.uiee) - 
    2 * adapw * (var.label - cov.label.uiee)
  
  res = list("label" = tvt, "var.label" = var.label, "uiee" = hvt, "var.uiee" = var.uiee, 
             "adaptive" = adaptive, "var.adaptive" = var.adaptive)
  
  return(res)
  
}



#################################################################################
# regression functions
cm0 = function(x) return(1 / 2)

cm1 = function(x, b)
{
  p = ncol(x)
  
  # b = rep(1, p)
  
  res = (x %*% b)
  
  return(res)
}



# cm2 = function(x, b)
# {
#   p = ncol(x)
# 
#   # b = rep(1, p)
#   delta = c(rep(0, p / 2), rep(1 / 3, p / 2))
# 
#   res = (x %*% b) * (1 - x %*% delta)
# 
#   return(res)
# }

cm2 = function(x, b)
{
  p = ncol(x)
  
  # b = rep(1, p)
  delta = c(rep(0, p / 2), rep(1 / 3, p / 2))
  
  res = (x %*% b) * (1 - x %*% delta)
  
  return(res)
}






cm3 = function(x, b)
{
  p = ncol(x)
  
  # b = rep(1, p)
  delta = c(rep(0, p / 2), rep(1 / 3, p / 2))
  
  res = (x %*% b) * (1 - x %*% delta) + sin(x) %*% rep(1, p)
  
  return(res)
}


cm4 = function(x, b)
{
  p = ncol(x)
  
  # b = rep(1, p)
  delta = c(rep(0, p / 2), rep(1 / 3, p / 2))
  
  res = (x %*% b) * (1 - x %*% delta) - exp( - (x %*% rep(1, p)) ) 
  
  return(res)
}



cm5 = function(x, b)
{
  p = ncol(x)
  
  # b = rep(1, p)
  delta = c(rep(0, p / 2), rep(1 / 3, p / 2))
  omega = rep(c(1 / 2, 0), p / 2)
  
  eta = rep(1, p)
  
  res = (x %*% b) * (1 - x %*% delta) - (x %*% omega) ^ 2 +
    log(abs(x %*% eta) + 1) * 2
  
  return(res)
}

#################################################################################

na.omit.list = function(y) { return(y[!vapply(y, function(x) any(is.na(x)), logical(1))]) }


