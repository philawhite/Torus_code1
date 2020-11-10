fit_NNGP = function(y,X_use,is_obs, times,hours, S, burn, mod_choice, thin=1, tune=100,lags = c(1,2,3,6,12,24*1:7) ){
  
  nt = length(times)
  ind_miss = which(is_obs == 0)
  n_miss = length(ind_miss)
  
  if(mod_choice == 0){
    neighs = lapply(1:nt,function(x){
      temp = x - 1; temp[temp >= 1] - 1
    })
  } else{
    neighs = lapply(1:nt,function(x){
      temp = x - lags; temp[temp >= 1] - 1
    })
  }
  
  
  U = inv_neigh(neighs, nt);
  
  if(mod_choice == 0){
    
    p = 1
    cand_var =  matrix(c(0.1),1,1)
    covs = numeric(p)
    covs[1] = 100 
    
    allow_prop = function(cand){
      ifelse( cand > 0 ,TRUE,FALSE)
    }
    
  } else if(mod_choice %in% c(7,8)){
    
    p = 1
    cand_var =  matrix(c(0.1),1,1)
    covs = numeric(p)
    covs[1] = 10
    
    allow_prop = function(cand){
      ifelse( cand > 0 ,TRUE,FALSE)
    }
    
  } else if(mod_choice == 1){
    
    p = 2
    cand_var = diag(c(1,10))/50
    covs = numeric(p)
    covs[1] = 10 ; covs[2] = 100
    
    allow_prop = function(q){
      ifelse( all(q[c(1,2)] > 0) & q[1] < 20 & q[2] < 200 ,TRUE,FALSE)
    }
    
  }  else if(mod_choice == 2){
    
    p = 4
    cand_var = diag(c(1,10,1,0.1))/50
    covs = numeric(p)
    covs[1] = 10 ; covs[2] = 100 ;covs[3] = 1 ;covs[4] = 1/2 
    
    allow_prop = function(q){
      ifelse( all(q[c(1,2,3,4)] > 0) & q[c(4)] <= 1,TRUE,FALSE)
    }
    
  }  else if(mod_choice == 3){
    
    p = 3
    cand_var = diag(c(1,10,0.1))/50
    covs = numeric(p)
    covs[1] = 10 ; covs[2] = 100 ;covs[3] = 1 
    
    allow_prop = function(q){
      ifelse( all(q[c(1,2,3)] > 0) & q[c(3)] <= 2,TRUE,FALSE)
    }
    
  } else if(mod_choice == 4){
    
    p = 5
    cand_var = diag(c(10,1,1,1,1))/50
    covs = numeric(p)
    covs[1] = 100 ; covs[2] = 1 ;covs[3] = 1 ;covs[4] = 0.1 ;covs[5] = 1
    
    allow_prop = function(q){
      ifelse( all(q[c(1,2,3,4,5)] > 0) & 
                (q[4] <= q[2]* q[3])  & 
                ((pi * (q[2] + q[3]) - pi^2 * q[4] ) > -1.0) ,TRUE,FALSE)
    }
    
  } else if(mod_choice == 5){
    
    p = 5
    cand_var = diag(c(1,10,1,1,1))/50
    covs = numeric(p)
    covs[1] = 10 ; covs[2] = 100 ;covs[3] = 1 ;covs[4] = 1/2; covs[5] = 1/2 
    
    allow_prop = function(q){
      ifelse( all(q[c(1,2,3,4,5)] > 0) & (q[5] <= 1) ,TRUE,FALSE)
    }
    
  } else if(mod_choice == 6){
    
    p = 4
    cand_var = diag(c(1,10,1,1))/50
    covs = numeric(p)
    covs[1] = 10 ; covs[2] = 100 ;covs[3] = 1 ;covs[4] = 1/2
    
    allow_prop = function(q){
      ifelse( all(q[c(1,2,3,4)] > 0) & (q[4] <= 1) ,TRUE,FALSE)
    }
    
  } else if(mod_choice == 9){
    
    p = 4
    cand_var = diag(c(0.01,0.01,0.1,1))/10
    covs = numeric(p)
    covs[1] =0.1 ; covs[2] = 0.1 ;covs[3] = 1 ;covs[4] = 100
    
    allow_prop = function(q){
      ifelse( all(q[c(1,2,3,4)] > 0) & sum(q[c(1,2)]) < 1 ,TRUE,FALSE)
    }
    
  }else if(mod_choice == 10){
    
    p = 4
    cand_var = diag(c(0.01,0.01,0.1,1))/10
    covs = numeric(p)
    covs[1] = 1 ; covs[2] = 1 ;covs[3] = 1 ;covs[4] = 100
    
    allow_prop = function(q){
      ifelse( all(q[c(1,2,3,4)] > 0) ,TRUE,FALSE)
    }
    
  }
  
  
  
  n_tot = S*thin + burn 
  
  a_t = 1
  b_t = 1
  a_s = 1
  b_s = 1
  
  
  ############ Impute, give initial values, and allot memory 
  
  sig2 = 1
  tau2 = 1
  
  X = cbind(1,X_use)
  xtx = t(X) %*% X
  q = ncol(X)
  
  Vbinv = solve(1e6 * diag(q))
  mb = rep(0,q)
  Vbinvmb = c(Vbinv %*% mb)
  
  bet = rep(0,q)
  bet[1] = mean(y,na.rm = TRUE)
  y[ind_miss] = bet[1]
  
  xb = X %*% bet
  w  = y - xb + rnorm(nt,0,0.5)
  
  res_sig2 = numeric(n_tot)
  res_tau2 = numeric(n_tot)
  res_bet = matrix(0,nrow =n_tot,ncol = q)
  res_w = matrix(0,nrow =n_tot,ncol = nt)
  res_covs = matrix(0,nrow =n_tot ,ncol = p)
  
  ############ Add an option for the exponential covariance
  
  time_difs = list_time_difs(times,neighs)
  day_difs = list_circ_difs_day(hours,neighs)
  week_difs = list_circ_difs_week(hours,neighs)
  
  cors = cor_list(day_difs,week_difs,time_difs,covs,mod_choice,nt) 
  Bs = Bmats( cors, nt)
  Fs = Fmats( cors, Bs, nt) 
  
  count = 0
  
  for(i in 2:n_tot){
    
    bet =  bet_update(Vbinv,xtx,Vbinvmb, w, tau2, y,X)
    xb = c(X %*% bet)
    
    tau2 = tau2_update(a_t , b_t, w,y,xb)
    sig2 = sig2_update(a_s , b_s, Bs,Fs, w, neighs)
    current_like = likelihood_compute( w, sig2, Fs , Bs, neighs)
    
    cand = mvrnorm(1,covs,cand_var) 
    
    if( allow_prop(cand) ){
      
      cov_out = cov_update(covs,cand,Bs,Fs, current_like,sig2,mod_choice, w,neighs,
                           count, day_difs,week_difs,time_difs)
      
      if(cov_out$count > count){
        
        covs = c(cov_out$covs)
        Bs = cov_out$B
        Fs = cov_out$F
        current_like = cov_out$like
        count = cov_out$count
        
      }
      
    }
    
    F_scale = Scale_F(Fs, sig2, nt)
    w =  c( w_update_all(w,is_obs,tau2, y ,  xb, neighs, U, Bs,  F_scale) )
    
    ############ Impute y
    
    y[ind_miss] = xb[ind_miss] + w[ind_miss] + rnorm(n_miss,0,sqrt(tau2))
    
    ########### Save output
    
    res_sig2[i] = sig2 
    res_tau2[i] = tau2
    res_bet[i,] = bet
    res_w[i,] = w 
    res_covs[i,] = covs 
    
    ##### Add tuner
    
    if(i %% tune == 0 & i <= burn){
      acc = count / tune
      count = 0
      
      if(acc > 0.02 & acc < 0.98){
        
        if(p > 1){
          temp_var = cov(res_covs[max(1,i - tune * 10 + 1):i,])
          rat = mean( diag(cand_var) / diag(temp_var) )
          cand_var = rat * temp_var 
        }
        
        if(acc < 0.1){
          cand_var = cand_var / 3
        }
        if(acc > 0.6){
          cand_var = cand_var * 2
        } 
      } else if(acc >= 0.98){
        cand_var = cand_var * 10
      } else{
        cand_var = cand_var / 10
      }
      
      cand_var = (cand_var + t(cand_var))/2 + 1e-6*diag(p)
    }
    
    ######## Iteration number
    if(i %% 10 == 0){
      cat("We are on: ", i,"\r") 
      flush.console() 
    }
    
  }
  
  ind_keep = which((1:n_tot > burn) & (1:n_tot %% thin == 0))
  
  out = vector(length = 5,mode = "list")
  out = list()
  # names(out) = c("sig2","tau2","w","mu","cov")
  
  out[["sig2"]] = res_sig2[ind_keep] ; rm(res_sig2)
  out[["tau2"]] =  res_tau2[ind_keep] ; rm(res_tau2)
  out[["w"]] = res_w[ind_keep,] ; rm(res_w)
  out[["bet"]] = res_bet[ind_keep,] ; rm(res_bet)
  out[["cov"]] = res_covs[ind_keep,] ; rm(res_covs)
  
  return(out)
  
}

