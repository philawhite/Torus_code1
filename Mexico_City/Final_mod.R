library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
#library(RcppSugar)
library(inline)
library(fields)
library(splancs)
library(coda)
library(xtable)
library(scoringRules)
library(MASS)
library(parallel)
#library(doMC)
library(foreach)

#RcppArmadillo.package.skeleton()

try(setwd("C:/Users/philaw/Box/Research/Porcu/Torus/code/ozone"),silent = TRUE)
rm(list=ls())
########################## Load Data

Rcpp::sourceCpp("nngp_circle.cpp")
source("fit_nngp.R")

in_function = function(x,val,alp=0.1){
  CI = quantile(x,c(alp/2,1 - alp/2))
  1*(val <= CI[2] & val >= CI[1])
}

O3_mat = read.csv("../../data/O3.csv")
PM10_mat = read.csv("../../data/PM10.csv")
RH_mat = read.csv("../../data/RH.csv")
TMP_mat = read.csv("../../data/TMP.csv")

times = 1:nrow(O3_mat)
#times = 1:(24*2*20) 

nt = length(times)
hours = O3_mat$HORA - 1

set.seed(1)
y = O3_mat$AJM
#y = PM10_mat$XAL

X_use = cbind(RH_mat$AJM,TMP_mat$AJM)



#ind_miss = which(is.na(y))
#ind_hold = sort(sample((1:(nt)),(nt)*0.2))
#ind_hold = sort(c(sapply(day_hold,function(x){ (x-1)*24 +(1:24)})))

#n_hold = length(ind_hold)
#y_holdout = y[ind_hold]
#y[ind_hold] = NA

S = 10000
burn = 10000
thin = 10


lags = c(1,2,3, 
         23,24,25,   
         47,48,49,  
         71,72,73,
         167,168,169,
         335,336,337,
         504, 1176)

is_obs =  1*(!is.na(y))
n_miss = sum(1 - is_obs)  
tune = 100

#comp_func = function(mmmm){

mod_choice = 9


mod = fit_NNGP(y = sqrt(y),X_use = X_use,is_obs =is_obs, times,hours, S = S, burn = burn, 
               mod_choice = mod_choice, thin=thin, tune=tune,lags = lags )

preds = matrix(0,ncol= nt,nrow = S)

for(i in 1:S){
  preds[i,] = mod$bet[i,1] + X_use %*% mod$bet[i,-1] + mod$w[i,] + rnorm(nt,0,sqrt(mod$tau2[i]))
}

preds = preds^2

#plot(times,y)
#lines(times,apply(preds,2,mean) ,col="black")
#lines(times,apply(preds,2,quantile,0.05) ,col="red",lty = 2)
#lines(times,apply(preds,2,quantile,0.95) ,col="red",lty = 2)

# pred = pred_NNGP(dat[,"time"],dat[,"lon"],dat[,"lat"],griddy[,"time"],griddy[,"lon"],
# griddy[,"lat"],neigh_grid,mod1$mu,mod1$sig2,mod1$tau2,mod1$w,mod1$cov,1)

# MAE = mean(abs( y_holdout - apply(preds[,ind_hold],2,mean) ))
# MSE = sqrt(mean(( y_holdout - apply(preds[,ind_hold],2,mean) )^2))
# crps = numeric(n_hold)
# covs = numeric(n_hold)
# 
# for(i in 1:n_hold){
#   covs[i] = in_function(preds[,ind_hold[i]], y_holdout[i] )
#   crps[i] = crps_sample(y_holdout[i] ,preds[,ind_hold[i]])
#   
# } 
# 
# mcrps = mean(crps)
# mcovs = mean(covs)
# 
# (out = c(mcrps,MAE,MSE,mcovs))
# 
# #}
# 
# 
# #out = mclapply( (0:11),function(x){ comp_func(x) },mc.cores = 4)
# #out = t(simplify2array(out) )
# #colnames(out) = c('CRPS',"MAE","MSE","90% Coverage")
# #xtable(out,digits = 4)
# 
 rm(list=setdiff(ls(), c("mod","preds")))

save.image("final_mode.RData")
