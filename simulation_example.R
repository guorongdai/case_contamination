# This is an example for the simulations of the case control project.
.libPaths("/home/guorongdai/RLibrary")

options(warn = - 1)

suppressMessages({
  if (!require(doParallel)) install.packages('doParallel', repos =
                                               'https://cran.revolutionanalytics.com/')
  library(doParallel)
  
  
  library(MASS)
})


source("ccs_functions.R")

registerDoParallel(detectCores())

j1 = 1
# j1 = 1 means the association model is logistic; j1 = 2 means the association is not logistic
n = 5000 # number of observations used in the study
nn = 50000 # raw sample size. 
# The sample used for the simulation is randomly selected from the raw data so that the numbers of
# cases and controls satisfy the requirement.
p = 12 # dimensionality of X_{-1}
cm = cm2 # regression function of the phenotyping model
nv = 200 # size of the validated set
m = 2 # Which component of theta is considered for inference?

############################################
# standard deviation of X
sigma = sqrt(1) # standard error of the normal distribution
############################################

############################################
# slope and intercept of the model generating raw data when the associaiton model is not logistic
beta = rep(1 / 2, p) 
beta0 = 1 / 2 
############################################

prop = 3 / 5 # proportion of controls in the study

b = rep(1, p) * 2 # slope of the linear term in the phenotyping model

############################################
# true value of the parameter when the association model is logistic
theta = rep(1, p) / 2
theta0 = 1 / 2
############################################


############################################
s = 500 # number of repetitions
out = numeric(4 * p + 12)

if(j1 == 1)
{
  
  ff = foreach(i = 1 : s, .combine = rbind) %dopar%
    {
      set.seed(512 * i)
      
      x.raw = matrix(rnorm(nn * p, 0, sigma), nrow = nn, ncol = p) # covariate
      
      pp = (x.raw %*% theta)[, 1] + theta0
      p0 = exp(cm(x.raw, b)[, 1]) # conditional probability of D = 0
      p1 = exp(pp) # conditional probability of D = 1
      p2 = rep(1, nn) # conditional probability of D = 2
      
      d.raw = numeric(nn)
      for(j in 1 : nn)
      {
        dd = rmultinom(1, 1, c(p0[j], p1[j], p2[j]))
        d.raw[j] = which(dd == 1) - 1
      }
      
      s.raw = (d.raw != 2)
      
      
      loc0 = which(s.raw == 0)
      loc1 = which(s.raw == 1)
      loc = c(loc0[1 : floor(n * prop)], loc1[1 : floor( n * (1-prop) )])
      
      n = length(loc)
      
      x = x.raw[loc, ]
      s = s.raw[loc]
      d = d.raw[loc]
      
      ind = which(s == 1) # locations of cases
      nca = length(ind) # size of the case pool
      
      
      r = rep(1, n) 
      r[ ind[-(1 : nv)] ] = 0
      # r = 0 means a case has not been tested; r = 1 means a case has been tested.
      
      validset = which(r == 1 & d != 0) # controls and real cases in the validation set
      
      ###################################
      # VEE
      res.valid = glm(s ~ x, family = binomial, subset = validset) 
      est.valid = coef(res.valid)[-1]
      sd.valid =  sqrt(diag(vcov(res.valid))[m])
      ci.valid = c(est.valid[m - 1] - 1.96 * sd.valid, est.valid[m - 1] + 1.96 * sd.valid)
      ###################################
      
      ###################################
      # NEE
      res.naive = glm(s ~ x, family = binomial) 
      est.naive = coef(res.naive)[- 1]
      sd.naive =  sqrt(diag(vcov(res.naive))[m])
      ci.naive = c(est.naive[m - 1] - 1.96 * sd.naive, est.naive[m - 1] + 1.96 * sd.naive)
      ###################################
      
      #############################
      # PIEE
      val = (r == 1) & (s == 1)
      inval = (r == 0) & (s == 1)
      d.val = d[val]
      x.val = x[val, ]
      x.inval = x[inval, ]
      
      piee.imp.mod = coef(glm(d.val ~ x.val, family = "binomial"))
      cond.exp.d = hf(addone(x.inval) %*% piee.imp.mod)[, 1]
      
      weights = rep(1, n)
      weights[inval] = cond.exp.d
      
      linear.set = which( r == 0 | (r == 1 & d != 0) ) 
      res.piee = suppressWarnings(glm(s ~ x, family = "binomial", weights = weights, subset = linear.set))
      est.piee = coef(res.piee)[- 1]
      sd.piee = sqrt( diag(piee.cov(coef(res.piee), piee.imp.mod, x, x, r, s, d))[m] )
      ci.piee = c(est.piee[m - 1] - 1.96 * sd.piee, est.piee[m - 1] + 1.96 * sd.piee)
      #############################
      
      #############################
      # UIEE
      sp = solver(x = x, w = x, d = d, s = s, r = r, sp = F, intercept1 = T, intercept2 = T, 
                  fi = T, eps = 1e-8, maxout = 20, intval = NULL) 
      est = sp $ theta # estimator
      hvar = sp $ covariance[m, m] # variance
      ci = c(est[m - 1] - 1.96 * sqrt(hvar), est[m - 1] + 1.96 * sqrt(hvar)) # 95% confidence interval
      #############################
      
      ###################################
      # bias
      out[1 : p] = est.valid - theta
      out[(p + 1) : (2 * p)] = est.naive - theta
      out[(2 * p + 1) : (3 * p)] = est.piee - theta
      out[(3 * p + 1) : (4 * p)] = est - theta
      ###################################
      
      ###################################
      # MSE
      out[4 * p + 1] = sum((est.valid - theta) ^ 2)
      out[4 * p + 2] = sum((est.naive - theta) ^ 2)
      out[4 * p + 3] = sum((est.piee - theta) ^ 2)
      out[4 * p + 4] = sum((est -theta) ^ 2)
      ###################################
      
      ###################################
      # CI
      out[4 * p + 5] = 1.96 * 2 * sd.valid
      out[4 * p + 6] = ! ( (theta[m] < ci.valid[1]) | (theta[m] > ci.valid[2]) )
      out[4 * p + 7] = 1.96 * 2 * sd.naive
      out[4 * p + 8] = ! ( (theta[m] < ci.naive[1]) | (theta[m] > ci.naive[2]) )
      out[4 * p + 9] = 1.96 * 2 * sd.piee
      out[4 * p + 10] = ! ( (theta[m] < ci.piee[1]) | (theta[m] > ci.piee[2]) )
      out[4 * p + 11] = 1.96 * 2 * sqrt(hvar)
      out[4 * p + 12] = ! ( (theta[m] < ci[1]) | (theta[m] > ci[2]) )
      ###################################
      
      out
    }
  
}


if(j1 == 2)
{
  
  binary = NULL # column number of binary covariates
  ############################################
  # standard deviation of X
  mu = numeric(p)
  sigma = diag(p) # standard error of the normal distribution
  ############################################
  
  ############################################
  beta = rep(1 / 2, p) # slope of the model generating raw data
  beta0 = 1 / 2 # intercept of the model generating raw data 
  ############################################
  
  prop = 3 / 5 # proportion of controls in the study
  
  ############################################
  # true value of the parameter approximated by a random sample
  nn0 = 200000 # total sample size
  n0 = nn0 / 2 # case-control sample size
  nca0 = 40000 # number of cases
  nc0 = n0 - nca0 # number of controls
  set.seed(77)
  x0.raw = mvrnorm(nn0, mu, sigma) # covariate
  x0.raw[, binary] = x0.raw[, binary] > 0 # make some binary covariates
  
  prob0 = hf( (x0.raw %*% beta)[, 1] + beta0 ) # conditional probability of S = 1
  s0.raw = rbinom(nn0, 1, prob0) # s = 0 means a control; s = 0 means a case.
  ind01 = which(s0.raw == 1) # locations of cases
  ind00 = which(s0.raw == 0) # locations of controls
  x0 = x0.raw[ c(ind01[1 : nca0], ind00[1 : nc0]), ] # case control sampling
  s0 = c( rep(1 , nca0) , rep(0, nc0) ) # case control sampling
  d0 = rep(2, n0)
  set.seed(7)
  d0[1 : nca0] = rbinom( nca0, 1, hf( cm( x0[1 : nca0, ], b) ) )
  # d = 0 means a false case; d = 1 means a real case; d = 2 means a control
  inclusion0 = which(d0 != 0) # observations which should be included in the study
  theta = glm(s0 ~ x0, family = binomial, subset = inclusion0) $ coefficients[-1]
  # true value of the parameter
  ############################################
  
  
  ############################################
  ff = foreach(i = 1 : s, .combine = rbind) %dopar%
    {
      set.seed(512 * i)
      
      x.raw = mvrnorm(nn, mu, sigma) # covariate
      x.raw[, binary] = x.raw[, binary] > 0 # make some binary covariates
      
      prob = hf( (x.raw %*% beta)[, 1] + beta0 ) # conditional probability of S = 1
      s.raw = rbinom(nn, 1, prob)  # s = 0 means a control; s = 1 means a case.
      
      loc0 = which(s.raw == 0)
      loc1 = which(s.raw == 1)
      loc = c(loc0[1 : floor(n * prop)], loc1[1 : floor( n * (1-prop) )])
      
      n = length(loc)
      
      x = x.raw[loc, ]
      s = s.raw[loc]
      
      ind = which(s == 1) # locations of cases
      nca = length(ind) # size of the case pool
      d = rep(2, n)
      set.seed(77 * i)
      d[ind] = rbinom( nca, 1, hf( cm(x[ind, ], b) ) )
      # d = 0 means a false case; d = 1 means a real case; d = 2 means a control
      
      inclusion = which(d != 0) # real cases and  controls
      
      
      r = rep(1, n) 
      r[ ind[-(1 : nv)] ] = 0
      # r = 0 means a case has not been tested; r = 1 means a case has been tested.
      
      validset = which(r == 1 & d != 0) # controls and real cases in the validation set
      
      ###################################
      # VEE
      res.valid = glm(s ~ x, family = binomial, subset = validset) 
      est.valid = coef(res.valid)[-1]
      sd.valid =  sqrt(diag(vcov(res.valid))[m])
      ci.valid = c(est.valid[m - 1] - 1.96 * sd.valid, est.valid[m - 1] + 1.96 * sd.valid)
      ###################################
      
      ###################################
      # NEE
      res.naive = glm(s ~ x, family = binomial) 
      est.naive = coef(res.naive)[- 1]
      sd.naive =  sqrt(diag(vcov(res.naive))[m])
      ci.naive = c(est.naive[m - 1] - 1.96 * sd.naive, est.naive[m - 1] + 1.96 * sd.naive)
      ###################################
      
      #############################
      # PIEE
      val = (r == 1) & (s == 1)
      inval = (r == 0) & (s == 1)
      d.val = d[val]
      x.val = x[val, ]
      x.inval = x[inval, ]
      
      piee.imp.mod = coef(glm(d.val ~ x.val, family = "binomial"))
      cond.exp.d = hf(addone(x.inval) %*% piee.imp.mod)[, 1]
      
      weights = rep(1, n)
      weights[inval] = cond.exp.d
      
      linear.set = which( r == 0 | (r == 1 & d != 0) ) 
      res.piee = suppressWarnings(glm(s ~ x, family = "binomial", weights = weights, subset = linear.set))
      est.piee = coef(res.piee)[- 1]
      sd.piee = sqrt( diag(piee.cov(coef(res.piee), piee.imp.mod, x, x, r, s, d))[m] )
      ci.piee = c(est.piee[m - 1] - 1.96 * sd.piee, est.piee[m - 1] + 1.96 * sd.piee)
      #############################
      
      #############################
      # UIEE
      sp = solver(x = x, w = x, d = d, s = s, r = r, sp = F, intercept1 = T, intercept2 = T, 
                  fi = T, eps = 1e-8, maxout = 20, intval = NULL) 
      est = sp $ theta # estimator
      hvar = sp $ covariance[m, m] # variance
      ci = c(est[m - 1] - 1.96 * sqrt(hvar), est[m - 1] + 1.96 * sqrt(hvar)) # 95% confidence interval
      #############################
      
      ###################################
      # bias
      out[1 : p] = est.valid - theta
      out[(p + 1) : (2 * p)] = est.naive - theta
      out[(2 * p + 1) : (3 * p)] = est.piee - theta
      out[(3 * p + 1) : (4 * p)] = est - theta
      ###################################
      
      ###################################
      # MSE
      out[4 * p + 1] = sum((est.valid - theta) ^ 2)
      out[4 * p + 2] = sum((est.naive - theta) ^ 2)
      out[4 * p + 3] = sum((est.piee - theta) ^ 2)
      out[4 * p + 4] = sum((est -theta) ^ 2)
      ###################################
      
      ###################################
      # CI
      out[4 * p + 5] = 1.96 * 2 * sd.valid
      out[4 * p + 6] = ! ( (theta[m] < ci.valid[1]) | (theta[m] > ci.valid[2]) )
      out[4 * p + 7] = 1.96 * 2 * sd.naive
      out[4 * p + 8] = ! ( (theta[m] < ci.naive[1]) | (theta[m] > ci.naive[2]) )
      out[4 * p + 9] = 1.96 * 2 * sd.piee
      out[4 * p + 10] = ! ( (theta[m] < ci.piee[1]) | (theta[m] > ci.piee[2]) )
      out[4 * p + 11] = 1.96 * 2 * sqrt(hvar)
      out[4 * p + 12] = ! ( (theta[m] < ci[1]) | (theta[m] > ci[2]) )
      ###################################
      
      out
    }
  
}

ff = na.omit(ff)

#############################
# bias
bias = ff[, 1 : (4 * p)]
bias = apply(bias, 2, mean)
vector.bias = numeric(4)
for(i in 1 : 4) vector.bias[i] = sqrt( sum(bias[(p * (i - 1) + 1) : (p * i)] ^ 2) ) / sqrt(p)
names(vector.bias) = c("VEE", "NEE", "PIEE", "UIEE")

print("####### Bias #####")
vector.bias
#############################

#############################
# MSE ratio
re = ff[, (4 * p + 1) : (4 * p + 4)]
re = apply(re, 2, mean)
re = re[1] / re
re = re[- 1]
names(re) = c("NEE", "PIEE", "UIEE")

print("####### MSE ratio #####")
re
#############################

#############################
# inference
# results of UIEE
inference = ff[, (4 * p + 5) : (4 * p + 12)]
inference = apply(inference, 2, mean)
names(inference) = c("CIL of VEE", "CR of VEE", "CIL of NEE", "CR of NEE",
                     "CIL of PIEE", "CR of PIEE", "CIL of UIEE", "CR of UIEE")

print("####### Inference #####")
inference
#############################
