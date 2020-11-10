fit_NNGP = function(y, NN, times,lon, lat, S, is_obs, burn, mod_choice, thin=1, tune=100){
  
  
  if(mod_choice == 5){  ## old 2
    
    p = 4
    cand_var = diag(c(0.01,0.01,0.1,0.1))/100
    covs = numeric(p)
    covs[1] = 5 ; covs[2] = 6 ;covs[3] = .005
    covs[4] = 0.9
    
    allow_prop = function(q){
      ifelse( all(q[c(1,2,3,4)] > 0) & q[c(4)] <= 1,TRUE,FALSE)
    }
    
  } else if(mod_choice == 6){ ## old 3
    
    p = 3
    cand_var = diag(c(0.1,0.1,0.1))/100
    covs = numeric(p)
    covs[1] = 3 ; covs[2] = 3 ;covs[3] = 1 
    
    allow_prop = function(q){
      ifelse( all(q[c(1,2,3)] > 0) & q[c(3)] <= 2,TRUE,FALSE)
    }
    
    
  } else if(mod_choice == 7){
    
    p = 5
    cand_var = diag(c(.1,.1,.1,.1,.1))/100
    covs = numeric(p)
    covs[1] = 5 ; covs[2] = 1 ;covs[3] = 1 ;covs[4] = 0.1 
    covs[5] = 1
    
    allow_prop = function(q){
      ifelse( all(q[c(1,2,3,4,5)] > 0) & 
                (q[4] <= q[2]* q[3])  & 
                ((pi * (q[2] + q[3]) - pi^2 * q[4] ) > -1.0) ,TRUE,FALSE)
    }
    
  } else if(mod_choice == 8){
    
    p = 5
    cand_var = diag(c(.1,.1,.1,.1,.1))/100
    covs = numeric(p)
    covs[1] = 10 ; covs[2] = 100 ;covs[3] = 1 ;covs[4] = 1/2
    covs[5] = 1/2
    
    allow_prop = function(q){
      ifelse( all(q[c(1,2,3,4,5)] > 0) & (q[5] <= 1) ,TRUE,FALSE)
    }
    
  } else if(mod_choice == 9){
    
    p = 4
    cand_var = diag(c(0.1,0.1,0.1,.1))/100
    covs = numeric(p)
    covs[1] = 0.3 ; covs[2] = 0.1 ;covs[3] = 1 
    covs[4] = 20
    
    allow_prop = function(q){
      ifelse( all(q[c(1,2,3,4)] > 0) & sum(q[c(1,2)]) < 1 ,TRUE,FALSE)
    }
    
  }else if(mod_choice == 10){
    
    p = 4
    cand_var = diag(c(0.1,0.1,0.1,.1))/100
    covs = numeric(p)
    covs[1] = 3 ; covs[2] = 3 ;covs[3] = 2 
    covs[4] = 3
    
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
  
  sig2 = 50
  tau2 = 50
  
  #X = cbind(1,X_use)
  #xtx = t(X) %*% X
  #q = ncol(X)
  
  #Vbinv = solve(1e6 * diag(q))
  #mb = rep(0,q)
  #Vbinvmb = c(Vbinv %*% mb)
  
  #bet = rep(0,q)
  #bet[1] = mean(y,na.rm = TRUE)
  #y[ind_miss] = bet[1]
  
  #xb = X %*% bet
  w  = y + rnorm(n,0,0.5)
  
  res_sig2 = numeric(n_tot)
  res_tau2 = numeric(n_tot)
  #res_bet = matrix(0,nrow =n_tot,ncol = q)
  res_w = matrix(0,nrow =n_tot,ncol = n)
  res_covs = matrix(0,nrow =n_tot ,ncol = p)
  
  ############ Add an option for the exponential covariance
  neighs = apply(NN,1,function(x){ 
    x[which(!is.na(x))] - 1
  } )
  
  time_difs = list_time_difs(dat_ordered$time,neighs)
  circ_difs = list_circ_difs_year(dat_ordered$time,neighs)
  dists = list_dists(dat_ordered$lon,dat_ordered$lat,neighs)
  
  cors = cor_list(dists,circ_difs,time_difs,covs,mod_choice,n) 
  Bs = Bmats( cors, n)
  Fs = Fmats( cors, Bs, n) 
  
  count = 0
  scale_fac = 2.38^2 / p
  
  st = proc.time()
  
  for(i in 2:n_tot){
    
    #  bet =  bet_update(Vbinv,xtx,Vbinvmb, w, tau2, y,X)
    #  xb = c(X %*% bet)
    
    tau2 = tau2_update(a_t , b_t, w,y)
    sig2 = sig2_update(a_s , b_s, Bs,Fs, w, neighs)
    current_like = likelihood_compute( w, sig2, Fs , Bs, neighs)
    
    cand = mvrnorm(1,covs,cand_var) 
    
    if( allow_prop(cand) ){
      
      cov_out = cov_update(covs,cand,Bs,Fs, current_like,sig2,mod_choice, w,neighs,
                           count, dists,circ_difs,time_difs)
      
      if(cov_out$count > count){
        
        covs = c(cov_out$covs)
        Bs = cov_out$B
        Fs = cov_out$F
        current_like = cov_out$like
        count = cov_out$count
        
      }
      
    }
    
    F_scale = Scale_F(Fs, sig2, n)
    w =  c( w_update_all(w,is_obs,tau2, y , neighs, U, Bs,  F_scale) )
    
    ############ Impute y
    
    # y[ind_miss] = xb[ind_miss] + w[ind_miss] + rnorm(n_miss,0,sqrt(tau2))
    
    ########### Save output
    
    res_sig2[i] = sig2 
    res_tau2[i] = tau2
    res_w[i,] = w 
    res_covs[i,] = covs 
    
    ##### Add tuner
    
    if(i %% tune == 0){
      
      acc = count / tune
      count = 0
      
      
      if(i > burn){
        
        cand_var = cov(res_covs[(burn/2):i,]) * scale_fac
        cand_var = (cand_var + t(cand_var))/2 + 1e-24*diag(p)
        
      } else{
        
        if(acc < 0.98 & acc > 0.02){
          
          if(i > burn/2){
            
            if(p > 1){
              cand_var = cov(res_covs[(burn/2):i,])
            }
            
          } else{
            
            if(p > 1){
              cand_var = cov(res_covs[max(1,i - tune * 10 + 1):i,])      
            }
            
            
          }
          
          if(acc < 0.05){
            scale_fac = scale_fac / 1.5
          }
          if(acc > 0.7){
            scale_fac = scale_fac * 2
          } 
        } else if(acc >= 0.98){
          scale_fac = scale_fac * 5
        } else{
          scale_fac = scale_fac / 10
        }
        
        cand_var = cand_var * scale_fac
        
        cand_var = (cand_var + t(cand_var))/2 + 1e-24*diag(p)
        
      }
      
    }  
    ######## Iteration number
    time_its <- (proc.time() - st)[3] / (i)
    time_used <- round((proc.time() - st)[3]/(60),digits=4)
    time_left <- round(time_its * (n_tot - i )/(60),digits=4)
    cat("\r", i, " of ", n_tot,"||| Time left: ",floor(time_left/60),
        " hours",time_left%%60," minutes") 
    flush.console()
    
  }
  
  ind_keep = 1:n_tot #which((1:n_tot > burn) & (1:n_tot %% thin == 0))
  
  out = vector(length = 4,mode = "list")
  out = list()
  # names(out) = c("sig2","tau2","w","mu","cov")
  
  out[["sig2"]] = res_sig2[ind_keep] ; rm(res_sig2)
  out[["tau2"]] =  res_tau2[ind_keep] ; rm(res_tau2)
  out[["w"]] = res_w[ind_keep,] ; rm(res_w)
  # out[["bet"]] = res_bet[ind_keep,] ; rm(res_bet)
  out[["cov"]] = res_covs[ind_keep,] ; rm(res_covs)
  
  return(out)
  
}

