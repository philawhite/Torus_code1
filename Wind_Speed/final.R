library(emulator)
library(parallel)
library(scoringRules)
library(fields)
library(MASS)

rm(list = ls())
try(setwd("C:/Users/philaw/Box/Research/Porcu/Torus/code/wind_backup/dic_comp"),silent = TRUE)

# Rcpp::sourceCpp("nngp_circle.cpp")

mod_choice = 4

pow = function(x,y){
  x^y
}

in_function = function(x,val,alp=0.1){
  CI = quantile(x,c(alp/2,1 - alp/2))
  1*(val <= CI[2] & val >= CI[1])
}

if(mod_choice == 1){
  
  p = 2
  cand_var = diag(c(1,10))/50
  covs = numeric(p)
  covs[1] = 1 ; covs[2] = 6
  
  allow_prop = function(q){
    ifelse( all(q[c(1,2)] > 0),TRUE,FALSE)
  }
  
  cov_func = function(theta,vartheta,time,pars){
    temp = exp(- vartheta / pars[1])
    res = exp(temp * cos(theta) - time / pars[2] - 1) * cos(temp * sin(theta))
    
    return(res)
  }
  
  prior_func = function(pars){
    sum(dgamma(pars[p + 1:2],0.01,0.01,log = TRUE)) + 
      dgamma(pars[1],1,1/pi,log = TRUE) + ### Seasonal distance c_s
      dgamma(pars[2],3,1/18,log = TRUE)  ### one year c_t
    
  }
  
}  else if(mod_choice == 2){
  
  p = 4
  cand_var = diag(c(1,10,1,0.1))/50
  covs = numeric(p)
  covs[1] = 1 ; covs[2] = 6 ;covs[3] = 1 ;covs[4] = 1/2 
  
  allow_prop = function(q){
    ifelse( all(q[c(1,2,3,4)] > 0) & q[c(4)] <= 1,TRUE,FALSE)
  }
  
  cov_func = function(theta,vartheta,time,pars){
    temp = exp(- vartheta / pars[1])
    res = pow((1 - pars[4]) / (1 - pars[4] * temp * cos(theta)),pars[3])
    
    return(res * exp(-time / pars[2]))
  }
  
  prior_func = function(pars){
    sum(dgamma(pars[p + 1:2],0.01,0.01,log = TRUE)) + 
      dgamma(pars[1],1,1/pi,log = TRUE) + ### Seasonal distance c_s
      dgamma(pars[2],3,1/18,log = TRUE) +  ### one year c_t
      dgamma(pars[3],1,1,log = TRUE)   ### \eta
    ### flat prior on epsilon
    
  }
  
  
}  else if(mod_choice == 3){
  
  p = 3
  cand_var = diag(c(1,10,0.1))/50
  covs = numeric(p)
  covs[1] = 1 ; covs[2] = 6 ;covs[3] = 1 
  
  allow_prop = function(q){
    ifelse( all(q[c(1,2,3)] > 0) & q[c(3)] <= 2,TRUE,FALSE)
  }
  
  cov_func = function(theta,vartheta,time,pars){
    temp1 = exp(- vartheta / pars[1])
    temp2 = pow(0.5,pars[3]) * (1 - temp1 * cos(theta))
    res = (1.0 - temp2) * exp(-time / pars[2])
    return(res)
  }
  
  prior_func = function(pars){
    sum(dgamma(pars[p + 1:2],0.01,0.01,log = TRUE)) + 
      dgamma(pars[1],1,1/pi,log = TRUE) + ### Seasonal distance c_s
      dgamma(pars[2],3,1/18,log = TRUE)   ### one year c_t
    ### flat prior on alpha
    
  }
  
  
} else if(mod_choice == 4){
  
  p = 5
  cand_var = diag(c(10,1,1,1,1))/50
  covs = numeric(p)
  covs[1] = 6 ; covs[2] = 1 ;covs[3] = 1 ;covs[4] = 0.1 ;covs[5] = 1
  
  allow_prop = function(q){
    ifelse( all(q[c(1,2,3,4,5)] > 0) & 
              (q[4] <= q[2]* q[3])  & 
              ((pi * (q[2] + q[3]) - pi^2 * q[4] ) > -1.0) ,TRUE,FALSE)
  }
  
  cov_func = function(theta,vartheta,time,pars){
    res = pow(1 + pars[2] * theta + pars[3] * vartheta - pars[4] * theta * vartheta, -pars[5])
    return(res * exp(-time / pars[1]))
  }
  
  prior_func = function(pars){
    sum(dgamma(pars[p + 1:2],0.01,0.01,log = TRUE)) + 
      dgamma(pars[1],3,1/18,log = TRUE) +  ### c_t
      dgamma(pars[5],1,1,log = TRUE)       ### \eta
    ### flat prior on all lambda
  }
  
} else if(mod_choice == 5){
  
  p = 5
  cand_var = diag(c(1,10,1,1,1))/50
  covs = numeric(p)
  covs[1] = 1 ; covs[2] = 6 ;covs[3] = 1 ;covs[4] = 1/2; covs[5] = 1/2 
  
  allow_prop = function(q){
    ifelse( all(q[c(1,2,3,4,5)] > 0) & (q[5] <= 1) ,TRUE,FALSE)
  }
  
  cov_func = function(theta,vartheta,time,pars){
    temp = 1 + vartheta / pars[1]
    res = 1 / pow( temp, pars[4] + pars[5]/2) * pow( 1 + (theta / pars[3]) / pow(temp,pars[5]), -1.0)
    return( res * exp(-time / pars[2]))
  }
  
  prior_func = function(pars){
    sum(dgamma(pars[p + 1:2],0.01,0.01,log = TRUE)) + 
      dgamma(pars[1],1,1/pi,log = TRUE) + ### Seasonal distance c_s
      dgamma(pars[2],3,1/18,log = TRUE) +  ### one year c_t
      dgamma(pars[3],1,1/pi,log = TRUE)   ### earth distance c_r
    dgamma(pars[4],1,1,log = TRUE)   ### \delta
    ## flat prior on b
  }
  
}  else if(mod_choice == 9){ ##### Example 6 in paper
  
  p = 4
  cand_var = diag(c(0.01,0.01,0.1,1))/10
  covs = numeric(p)
  covs[1] = 0.9 ; covs[2] = 0.05 ;covs[3] = 1 ;covs[4] = 6
  
  allow_prop = function(q){
    ifelse( all(q[c(1,2,3,4)] > 0) & sum(q[c(1,2)]) < 1 ,TRUE,FALSE)
  }
  
  cov_func = function(theta,vartheta,time,pars){
    temp1 = 1 - pars[1] - pars[2]
    temp2 = 1 - pars[1] * cos(theta) - pars[2]*cos(vartheta)
    res = pow(temp1,pars[3]) / pow(temp2, pars[3])
    return(res * exp(-time / pars[4]))
  }
  
  prior_func = function(pars){
    sum(dgamma(pars[p + 1:2],0.01,0.01,log = TRUE)) + 
      dgamma(pars[4],3,1/18,log = TRUE) +  ### one year c_t
      dgamma(pars[3],1,1,log = TRUE)   ### \eta
    ## flat prior on p_1 and p_2, subject to constraints
  }
  
} else if(mod_choice == 11){  ##### Separable Example in paper
  
  p = 3
  cand_var = diag(c(0.01,0.01,1))/10
  covs = numeric(p)
  covs[1] = 1 ; covs[2] = 1 ;covs[3] = 6
  
  allow_prop = function(q){
    ifelse( all(q[c(1,2,3)] > 0) ,TRUE,FALSE)
  }
  
  cov_func = function(theta,vartheta,time,pars){
    res = exp(- time/pars[3] - theta/pars[1] - vartheta/pars[2])
    return(res)
  }
  
  prior_func = function(pars){
    sum(dgamma(pars[p + 1:2],0.01,0.01,log = TRUE)) + 
      dgamma(pars[2],1,1/pi,log = TRUE) + ### Seasonal distance c_s
      dgamma(pars[3],3,1/18,log = TRUE) +  ### one year c_t
      dgamma(pars[1],1,1/pi,log = TRUE)   ### earth distance c_r
  }
  
}

prec_ou = function(difs, phi, sig2, n){
  
  prec = matrix(0,ncol = n,nrow = n)
  diags = numeric(n)
  
  if(length(unique(difs))==1){
    
    dif_u = unique(difs)
    diags[c(1,n)] = 1 / (1 - exp(-2 * phi *dif_u ))
    diags[2:(n-1)] = 1 / (1 - exp(-2 * phi *dif_u)) + exp(-2 * phi *dif_u) / 
      (1 - exp(-2 * phi * dif_u))
    
    diag(prec) = diags
    
    cross = - exp(-phi * dif_u) / (1 - exp(-2 * phi * dif_u))
    
    for(i in 1:(n-1)){
      prec[i,i+1] <- prec[i+1,i] <- cross
    }
    
  } else{
    
    diags[1] = 1 / (1 - exp(-2 * phi * difs[1] ))
    diags[2:(n-1)] = 1 / (1 - exp(-2 * phi * difs[1:(n-2)])) +
      exp(-2 * phi * difs[2:(n-1)]) / (1 - exp(-2 * phi * difs[2:(n-1)]))
    diags[n] = 1 / (1 - exp(-2 * phi * difs[n-1] ))
    
    diag(prec) = diags
    
    for(i in 1:(n-1)){
      prec[i+1,i] <- prec[i,i+1] <- - exp(-phi * difs[i]) / (1 - exp(-2 * phi * difs[i]))
    }
    
  }
  
  
  return(prec/sig2)
}


dmatnorm = function(X,mu,U_inv,V_inv,log_det_U,log_det_V){
  n = nrow(X)
  p = ncol(X)
  X_scale = sweep(X,2,mu,"-")
  left = V_inv %*% t(X_scale)
  right = U_inv %*% X_scale
  -n * p * log(2 * pi) / 2  - n/2 * log_det_V - p/2 * log_det_U - 1/2 * sum(diag(left %*% right))
}


my_dmvnorm = function(Y,sigma2,log_det,prec){
  n_local = length(Y) 
  -n_local/2 * log(2 * pi) - 
    1/2 * (log_det + n_local * log(sigma2)) - 
    quad.form(prec,Y) / (2 * sigma2)
}

my_dmvlnorm = function(logY,logmean,sigma2,log_det,prec){
  n_local = length(logY) 
  -n_local/2 * log(2 * pi) - 
    1/2 * (log_det + n_local * log(sigma2)) - 
    quad.form(prec,logY - logmean) / (2 * sigma2) - sum(logY)
}


likelihood = function(Y,log_det,prec){
  n_local = length(Y) 
  -n_local/2 * log(2 * pi) - 
    1/2 * (log_det) - 
    quad.form(prec,Y) / 2
}

dat_use = read.csv("wind_speed_use.csv")

n = nrow(dat_use)

set.seed(1)
y = dat_use$log_winds

X_use = cbind(dat_use$airt_max,dat_use$rh_max,dat_use$solar,
              sin(dat_use$hour * 2 * pi / 24),cos(dat_use$hour * 2 * pi / 24),
              sin(dat_use$hour * 2 * pi / 12),cos(dat_use$hour * 2 * pi / 12),
              sin(dat_use$hour * 2 * pi / 8),cos(dat_use$hour * 2 * pi / 8)
)
#out =  auto.arima(resid(lm(y~X_use)),d = 0)

#################################################
################ Hold out some days #############
#################################################

idx_all = sort(unique(c(which(is.na(y)) )))


n_hold = length(idx_all)

Y = y[-idx_all]
Y_hold = y[idx_all]
X = cbind(1,X_use[-idx_all,])
X_pred = cbind(1,X_use[idx_all,])


# with(dat_use[dat_use$month == 1,],
# points(proj_x,proj_y,cloud_cover_anom,legend.width=1, 
#        legend.shrink=1, horizontal=FALSE,xlab = "",ylab="",xaxt="n",yaxt = "n",
#        col = larry.colors(),nx = 30,ny = 30)
# )
# map('world',add=TRUE,col = "black",projecti on = "mollweide")
n = length(Y)
n_pred = length(Y_hold)

dd_angle = acos(cos(rdist(dat_use$windd[-idx_all] * 2 * pi / 360)))
dd_season = acos(cos(rdist(dat_use$time_numeric[-idx_all] * 2 * pi)))
dd = rdist(dat_use$time_numeric[-idx_all])

gd_angle = acos(cos(rdist(dat_use$windd[idx_all] * 2 * pi / 360,
                          dat_use$windd[-idx_all] * 2 * pi / 360)))
gd_season = acos(cos(rdist(dat_use$time_numeric[idx_all] * 2 * pi,
                           dat_use$time_numeric[-idx_all] * 2 * pi)))
gd = rdist(dat_use$time_numeric[idx_all],dat_use$time_numeric[-idx_all])

gg_angle = acos(cos(rdist(dat_use$windd[idx_all] * 2 * pi / 360)))
gg_season = acos(cos(rdist(dat_use$time_numeric[idx_all] * 2 * pi)))
gg = rdist(dat_use$time_numeric[idx_all])

# poly_terms1 = 8
# poly_terms2 = 6
# 
# n_poly = poly_terms1 * poly_terms2
# 
# gegenbauer_list1 = gegenbauer.polynomials(poly_terms1 - 1,1/2,normalized = TRUE)
# gegenbauer_list2 = gegenbauer.polynomials(poly_terms2 - 1,0,normalized = TRUE)


lm_temp = lm(Y ~  X[,-1])

q = ncol(X)


############################################################# 
############################################################# 
################## model 
############################################################# 
############################################################# 


reps = 2e4
burn = 5e4
tune = 100

b = array(0,c(reps,p)); b_now = covs
phi = rep(1/6,reps ); phi_now = 1/6
tau2 = rep(0.1,reps); tau2_now = 0.1
sig2 = rep(0.5,reps ) ; sig2_now = 0.5
beta_reg = matrix(0,reps,q) ; beta_reg_now =rep(0,q)
pred = matrix(0,reps,n_pred) 

#V_reg = array(0,c(reps,q,q)) ; V_reg_now = matrix(0,q,q)
like_save = numeric(reps)

beta_reg_now = coef(lm_temp)


pars_now = c(b_now,sig2_now,tau2_now)
pars_save = matrix(0,reps + burn,p + 2)
pars_save[1,] = pars_now

xb = X %*% beta_reg_now


# inv_R = prec_ou(rep(1,n_t),phi_beta_now,1,n_t)
# R = solve(inv_R)
# log_det_R = 2 * sum(log(diag(chol(R))))
# R_row_sum = apply(inv_R,1,sum)
# R_quad_form = sum(inv_R)


Sig_now = pars_now[p+1] * cov_func(dd_angle,dd_season,dd,pars_now[1:p]) +  
  pars_now[p+2] * diag(n)

chol_sig_now = chol(Sig_now)
log_sig_det = 2 * sum(log(diag(chol_sig_now)))
inv_sig = chol2inv(chol_sig_now)

XsigX = quad.form(inv_sig,X)

like_now = likelihood(Y - xb,log_sig_det,inv_sig)


scale_fac1 = 2.38^2 / (p+2)
cand_keep = 0.0001 * diag(p+2)
# cand_keep[1:p,1:p] = cand_var/ 10
# cand_keep[p+1,p+1] = 0.01
# cand_keep[p+2,p+2] = 0.01


cand_var1 = scale_fac1 * cand_keep


count1 = 0 ; chol1 = chol(cand_var1); inv1 = chol2inv(chol1); log_det1 = 2 * sum(log(diag(chol1)))

cand_var_phi = 0.1
count_phi = 0

st = proc.time()

for(i in 2:(reps + burn)){
  
  ######################### update Beta
  
  v_bet = solve(XsigX + diag(q)/100)
  m_bet = t(X) %*% c(inv_sig %*% Y)
  
  beta_reg_now = mvrnorm(1,v_bet %*% m_bet, v_bet)
  xb = X %*% beta_reg_now
  
  like_now = likelihood(Y - xb,log_sig_det,inv_sig)
  
  ######################### update all pars
  
  
  cand = c(pars_now + t(chol1) %*% rnorm(p+2))
  
  if(allow_prop(cand) & all(cand[-(1:p)] >0)){
    
    Sig_cand = cand[p+1] * cov_func(dd_angle,dd_season,dd,cand[1:p]) +  
      cand[p+2] * diag(n)
    
    chol_sig_cand= try(chol(Sig_cand),silent = TRUE)
    
    if(is.matrix(chol_sig_cand)){
      
      log_sig_det_cand = 2 * sum(log(diag(chol_sig_cand)))
      inv_sig_cand = chol2inv(chol_sig_cand)
      
      like_cand = likelihood(Y - xb,log_sig_det_cand,inv_sig_cand)
      
      prior_dif = prior_func(cand) - prior_func(pars_now)
      
      if(like_cand - like_now + prior_dif > log(runif(1))){
        
        pars_now = cand 
        
        Sig_now = Sig_cand
        log_sig_det= log_sig_det_cand
        inv_sig = inv_sig_cand
        
        like_now = like_cand
        XsigX = quad.form(inv_sig,X)
        count1 = count1 + 1
      }
      
    }
  }
  
  pars_save[i,] = pars_now
  
  
  
  if(i > burn){
    
    b[i - burn,] = pars_now[1:p]
    sig2[i - burn] = pars_now[p+1]
    tau2[i - burn] = pars_now[p+2]
    
    beta_reg[i - burn,] = beta_reg_now
    like_save[i - burn] = like_now
    
    
    #### Make Predictions
    Sig_gg = pars_now[p+1] * cov_func(gg_angle,gg_season,gg,pars_now[1:p]) + pars_now[p+2] * diag(n_pred) 
    Sig_gd = pars_now[p+1] * cov_func(gd_angle,gd_season,gd,pars_now[1:p])
    proj = Sig_gd %*% inv_sig
    cov_gg = Sig_gg - proj %*% t(Sig_gd)
    xb_pred = X_pred %*% beta_reg_now
    pred_mean = xb_pred + proj %*% (Y - xb)
    
    pred_temp = try(c(mvrnorm(1,pred_mean,cov_gg,tol = 1)),silent = TRUE)
    
    if(is.numeric(pred_temp)){
      pred[i - burn,] = pred_temp
    } else{
      pred[i - burn,] = NA
    }    
  }  
  
  
  if(i %% tune == 0){
    
    if(i < burn){
      
      acc1 = count1 / tune; count1 = 0
      
      if(acc1 > 0.05 & acc1 < 0.9){
        
        if(acc1 < 0.1){
          scale_fac1 = scale_fac1 / 1.5
        }
        
        if(acc1 > 0.6){
          scale_fac1 = scale_fac1 * 1.25
          
        }
        
      } else if(acc1 < 0.05){
        
        scale_fac1 = scale_fac1 / 3
        
      } else{
        
        scale_fac1 = scale_fac1 * 2
        
      }
      
      
      if(i > burn/10){
        
        cand_var1 = scale_fac1 * cov(unique(pars_save[1:i,]))
        cand_var1 = cand_var1 + diag(diag(cand_var1))/10
        
      } else{
        
        cand_var1 = scale_fac1 * cand_keep
        
      }
      
      
    } else{
      
      cand_var1 = scale_fac1 * cov( unique(pars_save[(burn/2):i,]) )
      cand_var1 = cand_var1 + diag(diag(cand_var1))/10
    }
    
    chol1 = chol(cand_var1); inv1 = chol2inv(chol1); log_det1 = 2 * sum(log(diag(chol1)))
    
    # time_its <- (proc.time() - st)[3] / (i)
    # time_used <- round((proc.time() - st)[3]/(60),digits=4)
    # time_left <- round(time_its * (reps + burn - i )/(60),digits=4)
    # cat("\r", i, " of ", reps + burn,"||| Time left: ",floor(time_left/60),
    #     " hours",time_left%%60," minutes") 
    # flush.console()
    # cat("\n,   acceptance rate:", acc1,"scale=",scale_fac1)
    
  }
  
  chol1 = chol(cand_var1); inv1 = chol2inv(chol1); log_det1 = 2 * sum(log(diag(chol1)))
  
  time_its <- (proc.time() - st)[3] / (i)
  time_used <- round((proc.time() - st)[3]/(60),digits=4)
  time_left <- round(time_its * (reps + burn - i )/(60),digits=4)
  cat("\r", i, " of ", reps + burn,"||| Time left: ",floor(time_left/60),
      " hours",time_left%%60," minutes")
  flush.console()
  cat("\n,   count1:", c(count1), "like_now= ",like_now,"scale_fac= ",scale_fac1 )
  
  
}

# devs = -2 * like_save
# 
# d_bar = mean(devs)
# 
# p_d = var(devs)/2
# 
# DIC = d_bar + p_d
# 
# (out = c(DIC,d_bar,p_d))

# 
# b[i - burn,] = pars_now[1:p]
# sig2[i - burn] = pars_now[p+1]
# tau2[i - burn] = pars_now[p+2]
# 
# beta_reg[i - burn,] = beta_reg_now
# like_save[i - burn] = like_now

rm(list=setdiff(ls(), c("like_save","b","sig2","tau2","beta_reg") ))

save.image("final.RData")
