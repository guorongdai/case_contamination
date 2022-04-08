rm(list = ls())

suppressMessages({
  
  if (!require(ggplot2)) install.packages('ggplot2', repos =
                                            'https://cran.revolutionanalytics.com/')
  
  if (!require(gridExtra)) install.packages('gridExtra', repos =
                                              'https://cran.revolutionanalytics.com/')
  
  
  library(ggplot2)
  library(gridExtra)
  
})

setwd("~/Dropbox/CCS/Sepsis_data")
source("ccs_functions.R")

dat = read.csv("sepsis_data.csv", header = T)
used.var = c("urineoutput", "lactate_min", "bun_mean", "sysbp_min", "age", 
             "sodium_max", "aniongap_max", "creatinine_min", "spo2_mean", 
             "is_male", "metastatic_cancer")
N = nrow(dat)

x = as.matrix(dat[, 2 : 12])
s = dat $ is_in_the_caes_pool
d = dat $ status
d[is.na(d)] = 1
r = dat $ is_validated



##########################################################
# compute various estimators
#############################
# UIEE
res.uiee = solver(x, x, d, s, r)
theta.uiee = res.uiee $ theta
var.uiee = diag(res.uiee $ covariance)[2 : 12]
#############################

#############################
# PIEE
val = (r == 1) & (s == 1)
inval = (r == 0) & (s == 1)
d.val = d[val]
x.val = x[val, ]
x.inval = x[inval, ]

piee.imp.mod = coef(glm(d.val ~ x.val, family = "binomial"))
cond.exp.d = hf(addone(x.inval) %*% piee.imp.mod)[, 1]

weights = rep(1, N)
weights[inval] = cond.exp.d

res.piee = suppressWarnings(glm(s ~ x, family = "binomial", weights = weights))
theta.piee = coef(res.piee)[-1]
var.piee = diag(piee.cov(coef(res.piee), piee.imp.mod, x, x, r, s, d))[2 : 12]
#############################

#############################
# naive
res.naive = glm(s ~ x, family = "binomial")
theta.naive = coef(res.naive)[-1]
var.naive = diag(vcov(res.naive))[-1]
#############################

#############################
# validation-only
res.valid = glm(s ~ x, family = "binomial", subset = ((r == 1) & (d != 0)))
theta.valid = coef(res.valid)[-1]
var.valid = diag(vcov(res.valid))[-1]
#############################
##########################################################

##########################################################
# generate the plot
p = ncol(x)
theta = c(theta.uiee, theta.piee, theta.naive, theta.valid)
sd = sqrt(c(var.uiee, var.piee, var.naive, var.valid))
CI_upper = theta + 1.96 * sd
CI_lower = theta - 1.96 * sd
method = as.factor(rep(c("UIEE", "PIEE", "NEE", "VEE"), each = p))
coef = as.factor( rep(used.var, 4) )

res = data.frame("theta" = theta, "sd" = sd, "CI_upper" = CI_upper,
                 "CI_lower" = CI_lower, "method" = method, "coef" = coef)

ggplot(res,aes(x=coef,y=theta,ymin=CI_lower,ymax=CI_upper,col=factor(method)))+
  geom_linerange(size=1,position=position_dodge(width=0.5))+coord_flip()+
  geom_pointrange(size=0.5,position=position_dodge(width=0.5)) +
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
  labs(y = NULL, x="Coefficient",col="Method",title=NULL)
##########################################################
