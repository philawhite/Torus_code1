#include <stdlib.h>            // malloc
#include <stdio.h>             // printf
#include <math.h>              // fabs, sqrt, etc.
#include <Rmath.h>              // fabs, sqrt, etc.
#include <time.h>              // time
#include <cmath>
#include <unistd.h>            // getpid
#include <string>
#include <RcppArmadillo.h>
//#include <Rcpp.h>
using namespace Rcpp;
using namespace R;
using namespace arma;  // use the Armadillo library for matrix computations

// [[Rcpp::export]]
int c_which(IntegerVector x, int look_for){
  int res;
  bool not_found = true;
  int i = 0;
  while(not_found){
    if(x[i] == look_for){
      res = i;
      not_found = false;
    }
    i = i + 1;
  }
  return res;
}

// [[Rcpp::export]]
double likelihood_compute(vec w, double sig2, List Fs , List Bs, List neigh){
  
  int n = w.size();
  double Ftemp = Fs[0];
  double new_like = dnorm(w[0],0,sqrt(sig2 * Ftemp),true);
  
  //  Rcpp::Rcout << "new like = " << new_like << std::endl;
  for(int i=1 ; i < n; i++){
    uvec t_neigh = neigh[i];
    vec B_temp = Bs[i];
    vec Bw =  B_temp.t() * w.elem(t_neigh) ;
    mat Ftemp = Fs[i];
    new_like += dnorm(w[i],Bw[0],sqrt(sig2 * Ftemp(0,0)),true);
    // Rcpp::Rcout << "new like = " << new_like << "||| we are on" << i << std::endl;
    
  }
  
  return new_like;
  
}


// double euc_sphere_dist(double long1,double lat1,double long2, double lat2){
//   double long1s = long1 * M_PI/180.0;
//   double long2s = long2 * M_PI/180.0;
//   double lat1s = lat1 * M_PI/180.0;
//   double lat2s = lat2 * M_PI/180.0;
//   
//   double x1 = cos(lat1s) * cos(long1s);
//   double y1 = cos(lat1s) * sin(long1s);
//   double z1 = sin(lat1s);
//   double x2 = cos(lat2s) * cos(long2s);
//   double y2 = cos(lat2s) * sin(long2s);
//   double z2 = sin(lat2s);
//   double out = sqrt(pow(x1 - x2,2.0) + pow(y1 - y2,2.0) + pow(z1 - z2,2.0) );
//   
//   return out;
// }

vec my_mvrnorm(vec mu, mat sigma){
  int ncols = sigma.n_cols;
  vec Y = rnorm(ncols,0,1);
  return mu + chol(sigma) * Y;
}

NumericMatrix euc_dist(NumericMatrix x){
  int n = x.nrow();
  double d;
  NumericMatrix out(n,n);
  
  for (int i = 0; i < n - 1; i++){
    NumericVector v1 = x.row(i);
    for (int j = i + 1; j < n ; j ++){
      d = sqrt(sum(pow(v1-x.row(j), 2.0)));
      out(j,i)=d;
      out(i,j)=d;
    }
  }
  
  return out;
}

NumericMatrix circ_dist_day(NumericVector x){
  int n = x.size();
  NumericMatrix out(n,n);
  
  for (int i = 0; i < n - 1; i++){
    double long1s = x[i] * M_PI/12.0;
    
    for (int j = i + 1; j < n ; j ++){
      
      double long2s = x[j] * M_PI/12.0;
      double dlong = fabs(long1s - long2s) ;
      double lon_use;
      
      if(dlong > M_PI){
        lon_use = 2*M_PI - dlong;
      } else{
        lon_use = dlong;
      }
      
      out(j,i) = lon_use;
      out(i,j) = lon_use;
    }
  }
  
  return out;
}


NumericMatrix circ_dist_week(NumericVector x){
  int n = x.size();
  NumericMatrix out(n,n);
  
  for (int i = 0; i < n - 1; i++){
    double long1s = x[i] * M_PI/84.0;
    
    for (int j = i + 1; j < n ; j ++){
      
      double long2s = x[j] * M_PI/84.0;
      double dlong = fabs(long1s - long2s) ;
      double lon_use;
      
      if(dlong > M_PI){
        lon_use = 2*M_PI - dlong;
      } else{
        lon_use = dlong;
      }
      
      out(j,i) = lon_use;
      out(i,j) = lon_use;
    }
  }
  
  return out;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List inv_neigh(List neighs,int n_all){
  List U(n_all);
  
  for(int i=0 ; i <  n_all; i++){
    IntegerVector inv_neigh;
    for(int j=1; j < n_all; j++){
      if(j != i){
        IntegerVector temp = neighs[j];
        bool b = any( temp == i ).is_true();
        if(b){
          inv_neigh.push_back(j);
        } 
      }
    }
    U[i] = inv_neigh;
  }
  
  return U;
}


////////////////////////////////////////////////////////////////////////
////////////////////// Distances /////////////////////////////////////
////////////////////////////////////////////////////////////////////////


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List list_time_difs(vec times,List neighs){
  int n = times.size();
  List dists(n);
  
  NumericMatrix first_dist(1,1);
  dists[0] = first_dist;
  
  for(int i=1 ; i < n ; i++){
    IntegerVector ind_neigh = neighs[i];
    int mat_size = ind_neigh.size() + 1;
    NumericMatrix times_neigh(mat_size,1); 
    times_neigh(0,0) = times[i];
    
    for(int j=1; j < mat_size; j++){
      times_neigh(j,0) = times[ind_neigh[j-1]]; 
    }
    
    NumericMatrix dist_mat(mat_size,mat_size);
    dist_mat = euc_dist(times_neigh);
    dists[i] = dist_mat;
  }
  return dists;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List list_circ_difs_day(vec hours,List neighs){
  int n = hours.size();
  List dists(n);
  
  NumericMatrix first_dist(1,1);
  dists[0] = first_dist;
  
  for(int i=1 ; i < n ; i++){
    
    IntegerVector ind_neigh = neighs[i];
    int mat_size = ind_neigh.size() + 1;
    NumericVector times_neigh(mat_size); 
    times_neigh[0] = hours[i];
    
    for(int j=1; j < mat_size; j++){
      times_neigh[j] = hours[ind_neigh[j-1]]; 
    }
    
    NumericMatrix dist_mat(mat_size,mat_size);
    dist_mat = circ_dist_day(times_neigh);
    dists[i] = dist_mat;
  }
  return dists;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List list_circ_difs_week(vec hours,List neighs){
  int n = hours.size();
  List dists(n);
  
  NumericMatrix first_dist(1,1);
  dists[0] = first_dist;
  
  for(int i=1 ; i < n ; i++){
    
    IntegerVector ind_neigh = neighs[i];
    int mat_size = ind_neigh.size() + 1;
    NumericVector times_neigh(mat_size); 
    times_neigh[0] = hours[i];
    
    for(int j=1; j < mat_size; j++){
      times_neigh[j] = hours[ind_neigh[j-1]]; 
    }
    
    NumericMatrix dist_mat(mat_size,mat_size);
    dist_mat = circ_dist_week(times_neigh);
    dists[i] = dist_mat;
  }
  return dists;
}

// List list_grid_time_difs(vec times,vec grid_times,List neighs){
//   int n = grid_times.size();
//   List dists(n);
//   for(int i=0 ; i < n ; i++){
//     IntegerVector ind_neigh = neighs[i];
//     int mat_size = ind_neigh.size() + 1;
//     NumericMatrix times_neigh(mat_size,1); 
//     times_neigh(0,0) = grid_times[i];
    
//     for(int j=1; j < mat_size; j++){
//       times_neigh(j,0) = times[ind_neigh[j-1]]; 
//     }
//     NumericMatrix dist_mat(mat_size,mat_size);
//     dist_mat = euc_dist(times_neigh);
//     dists[i] = dist_mat;
//   }
//   return dists;
// }

// List list_grid_circ_difs(vec hours,vec grid_hours,List neighs){
//   int n = grid_hours.size();
//   List dists(n);
//   for(int i=0 ; i < n ; i++){
//     IntegerVector ind_neigh = neighs[i];
//     int mat_size = ind_neigh.size() + 1;
//     NumericVector times_neigh(mat_size); 
//     times_neigh[0] = grid_hours[i];
    
//     for(int j=1; j < mat_size; j++){
//       times_neigh[j] = hours[ind_neigh[j-1]]; 
//     }
//     NumericMatrix dist_mat(mat_size,mat_size);
//     dist_mat = circ_dist(times_neigh);
//     dists[i] = dist_mat;
//   }
//   return dists;
// }

////////////////////////////////////////////////////////////////////////
////////////////////// Covariances /////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
double exp_decay(double time,double ct){
  double res = exp(-time/ct);
  return res;
}

// [[Rcpp::export]]
double exp_circle_day(double time,double ch){
  double res = exp(-time/ch);
  return res;
}

// [[Rcpp::export]]
double exp_circle_week(double time,double cw){
  double res = exp(-time/cw);
  return res;
}

// [[Rcpp::export]]
double example1(double theta,double vartheta,double time,
double cw, double ct){
  double temp = exp(- vartheta / cw);
  double res = exp(temp * cos(theta) - time / ct - 1) * cos(temp * sin(theta));

  return res;
}

// [[Rcpp::export]]
double example2(double theta,double vartheta,double time, 
double cw, double ct, double tau, double eps){
  double temp = exp(- vartheta / cw);

  double res = pow((1 - eps) / (1 - eps * temp * cos(theta)),tau);
  return res * exp(-time / ct);
}

// [[Rcpp::export]]
double example3(double theta,double vartheta,double time,
double cw, double ct, double alp){
  double temp1 = exp(- vartheta / cw);
  double temp2 = pow(0.5,alp) * (1 - temp1 * cos(theta));
  double res = (1.0 - temp2) * exp(-time / ct);
  return res;
}

// [[Rcpp::export]]
double example4(double theta,double vartheta,double time,
double ct, double lam1, double lam2, double lam3, double tau){
  double res = pow(1 + lam1 * theta + lam2 * vartheta - lam3 * theta * vartheta, -tau);
  return res * exp(-time / ct);
}

// [[Rcpp::export]]
double example5(double theta,double vartheta,double time,
double cw, double ct, double ch, double del, double bet){
  double temp = 1 + vartheta / cw;
  double res = 1 / pow( temp, del + bet/2) * pow( 1 + (theta / ch) / pow(temp,bet), -1.0);
  return res * exp(-time / ct);
}

// [[Rcpp::export]]
double example6(double theta,double vartheta,double time,
double cw, double ct, double ch, double bet){
  double temp = 1 + vartheta / cw;
  double res = 1 / pow( temp, 2.0 * bet) * exp(-(theta / ch) * pow(temp,bet));
  return res * exp(-time / ct);
}

// [[Rcpp::export]]
double example9(double theta,double vartheta,double time,
double p1, double p2, double tau, double ct){
  double temp1 = 1 - p1 - p2;
  double temp2 = 1 - p1 * cos(theta) - p2*cos(vartheta);
  double res = pow(temp1,tau) / pow(temp2, tau);
  return res * exp(-time / ct);
}

// [[Rcpp::export]]
double example10(double theta,double vartheta,double time,
double lam1, double lam2, double lam3, double ct){
  double temp1 = (cos(theta) -1) ;
  double temp2 = (cos(vartheta) - 1);
  double res = lam1 * temp1 + lam2 * temp2 + lam3 * temp1 * temp2;
  return exp(res - time / ct);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List cor_list(List circ_day, List circ_week, List times,vec pars,int mod_choice,int n){
  List cov(n); 
  double one_cov = 1;
  
  for(int i=0 ; i < n; i++){
    mat theta =as<arma::mat>(circ_day[i]);
    mat vartheta =as<arma::mat>(circ_week[i]);
    mat temp_difs =as<arma::mat>(times[i]);
    
    int n_mat = theta.n_rows;
    mat cov_dist(n_mat,n_mat);
    
    for(int j = 0; j < n_mat ; j++){
      for(int k = j ; k < n_mat ; k++){
        
        if(mod_choice == 0){  // Simple Model
          one_cov =  exp_decay(temp_difs(j,k),
          pars[0]);
        }
        if(mod_choice == 7){  // Simple Model
          one_cov =  exp_circle_day(theta(j,k),
          pars[0]);
        }
        if(mod_choice == 8){  // Simple Model
          one_cov =  exp_circle_week(vartheta(j,k),
          pars[0]);
        }
        if(mod_choice == 1){  // Example 1
          one_cov =  example1(theta(j,k),vartheta(j,k),temp_difs(j,k),
          pars[0], pars[1]);
        }
        if(mod_choice == 2){ // Example 2
          one_cov = example2(theta(j,k),vartheta(j,k),temp_difs(j,k),
          pars[0],pars[1],pars[2],pars[3]);
        }
        if(mod_choice == 3){ // Example 2
          one_cov = example3( theta(j,k),vartheta(j,k),temp_difs(j,k),
          pars[0],pars[1],pars[2]);
        }
        if(mod_choice == 4){  // Example 2
          one_cov = example4(theta(j,k),vartheta(j,k),temp_difs(j,k),
          pars[0],pars[1],pars[2],pars[3],pars[4]);
        }
        if(mod_choice == 5){ // Example 2
          one_cov = example5(theta(j,k),vartheta(j,k),temp_difs(j,k),
          pars[0],pars[1],pars[2],pars[3],pars[4]);
        }
        if(mod_choice == 6){ // Example 2
          one_cov = example6(theta(j,k),vartheta(j,k),temp_difs(j,k),
          pars[0],pars[1],pars[2],pars[3]);
        }
        if(mod_choice == 9){ // Example 2
          one_cov = example9(theta(j,k),vartheta(j,k),temp_difs(j,k),
          pars[0],pars[1],pars[2],pars[3]);
        }

        if(mod_choice == 10){ // Example 2
          one_cov = example10(theta(j,k),vartheta(j,k),temp_difs(j,k),
          pars[0],pars[1],pars[2],pars[3]);
        }
        
        cov_dist(j,k)= one_cov;
        cov_dist(k,j)= one_cov;
      }
    }
    
    cov[i] = cov_dist;
  }
  
  return cov; 
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List Bmats(List covs,int n){
  List Bmat(n);
  Bmat[0] = 0;
  for(int i=1; i < n ; i++){
    // Rcpp::Rcout << "i = " << i << std::endl;
    mat temp = covs[i];
    int mat_size = temp.n_rows;
    mat C_s = temp.submat( span(0, 0) , span(1, mat_size-1) ) ;
    mat C_N = temp.submat( span(1, mat_size-1) , span(1, mat_size-1) ) +  1e-6 * eye<mat>(mat_size-1,mat_size-1) ;
    mat res = C_s * C_N.i();
    Bmat[i] = res;
  }
  return Bmat;
}

// [[Rcpp::depends(RcppArmadillo)]]
List Bmats_g(List covs,int n){
  List Bmat(n);
  for(int i=0; i < n ; i++){
    mat temp = covs[i];
    int mat_size = temp.n_rows;
    mat C_s = temp.submat( span(0, 0) , span(1, mat_size-1) ) ;
    mat C_N = temp.submat( span(1, mat_size-1) , span(1, mat_size-1) ) +  1e-12 * eye<mat>(mat_size-1,mat_size-1);
    mat res = C_s * C_N.i();
    Bmat[i] = res;
  }
  return Bmat;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List Fmats(List covs,List Bmats,int n){
  List Fmat(n);
  Fmat[0] = 1;
  for(int i=1; i < n ; i++){
    mat temp = covs[i];
    mat B = Bmats[i];
    int mat_size = temp.n_rows;
    double C_s = temp(0,0) ;
    mat C_sn = temp.submat( span(1, mat_size-1) , span(0,0) ) ;
    mat res = C_s - B * C_sn;
    Fmat[i] = res;
  }
  return Fmat;
}

// [[Rcpp::depends(RcppArmadillo)]]
List Fmats_g(List covs,List Bmats,int n){
  List Fmat(n);
  for(int i=0; i < n ; i++){
    mat temp = covs[i];
    mat B = Bmats[i];
    int mat_size = temp.n_rows;
    double C_s = temp(0,0) ;
    mat C_sn = temp.submat( span(1, mat_size-1) , span(0,0) ) ;
    mat res = C_s - B * C_sn;
    Fmat[i] = res;
  }
  return Fmat;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List Scale_F(List F,double sig2,int n){
  List SF(n);
  SF[0] = sig2;
  for(int i=1; i < n ; i++){
    mat Ftemp = F[i];
    SF[i] = sig2*Ftemp;
  }
  return SF;
}


// [[Rcpp::depends(RcppArmadillo)]]
vec a_sum(List neigh,IntegerVector U,vec w,List B, List F,int obs_ind){ 
  NumericVector mV(2);
  int n_U = U.size();
  double m_term = 0 ;
  double V_term = 0 ;
  
  for(int i=0 ; i< n_U; i++ ){   //// t in U(s)
    double w_temp = w[U[i]];
    
    IntegerVector t_neigh = neigh[U[i]];
    int n_t = t_neigh.size();
    int obs_in_neigh = c_which(t_neigh,obs_ind);
    vec Btemp = B[U[i]]; /// B_{t}
    mat Ftemp = F[U[i]]; /// F_{t}
    
    
    for(int j=0; j < n_t; j++){  //// s in N(t)
      if(j != obs_in_neigh){
        w_temp += - Btemp[j] * w[t_neigh[j]] ;
      }
    }
    m_term += Btemp[obs_in_neigh] * w_temp/ Ftemp(0,0);
    V_term += pow(Btemp[obs_in_neigh],2.0)/ Ftemp(0,0);
  }
  
  mV[0] = m_term;
  mV[1] = V_term;
  
  return(mV);
}

// [[Rcpp::depends(RcppArmadillo)]]
double w_update(int obs_ind,bool is_observed,vec &w,double tau2,vec y , 
vec xb,List neigh,IntegerVector U,List B, List F){ 
  
  int n = y.size();
  uvec neigh_ind = neigh[obs_ind];
  vec U_terms = a_sum(neigh, U, w, B,  F, obs_ind) ;
  
  if(obs_ind == (n-1)){
    U_terms[0] = 0;
    U_terms[1] = 0;
  }
  
  mat F_w = F[obs_ind];
  vec B_w = B[obs_ind];
  vec w_N = w.elem(neigh_ind);
  
  //  double V_w = 1/(1/F_w(0,0) + U_terms[1]  + 1/tau2) ;
  //  mat m_w =  B_w.t() * w_N / F_w(0,0) + U_terms[0] + (y[obs_ind] - xb[obs_ind])/tau2;
  
  double V_w; mat m_w(1,1);
  
  if(is_observed){
    V_w = 1 / ( 1/F_w(0,0) + U_terms[1]  + 1/tau2 );
    m_w = B_w.t() * w_N / F_w(0,0) + U_terms[0] + (y[obs_ind] - xb[obs_ind])/tau2;
  } else{
    V_w = 1 / ( 1/F_w(0,0) + U_terms[1] ) ;
    m_w =B_w.t() * w_N / F_w(0,0) + U_terms[0] ;
  }
  
  w[obs_ind] = rnorm(1,V_w*m_w(0,0),sqrt(V_w))[0];
  return w[obs_ind];
  
}

// [[Rcpp::depends(RcppArmadillo)]]
double w_update0(int obs_ind,bool is_observed,vec &w,double tau2,vec y , 
vec xb,List neigh,IntegerVector U,List B, List F){ 
  
  uvec neigh_ind = neigh[obs_ind];
  
  vec U_terms = a_sum(neigh, U, w, B,  F, obs_ind) ;
  
  double F_w = as<double>(F[obs_ind]);
  double V_w; double m_w;
  
  if(is_observed){
    V_w = 1 /( 1/F_w + U_terms[1] + 1/tau2) ;
    m_w = U_terms[0] + (y[obs_ind] - xb[obs_ind])/tau2;
  } else{
    V_w = 1 / ( 1/F_w + U_terms[1] ) ;
    m_w = U_terms[0] ;
  }
  
  w[obs_ind] = rnorm(1,V_w*m_w,sqrt(V_w))[0];
  return w[obs_ind];
  
}

// [[Rcpp::export]]
vec w_update_all(vec w,IntegerVector is_obs, double tau2,vec y , vec xb,List neigh,
                 List U,List B, List F_scale){
  
  int n = w.size();
  vec out(n);
  
  out[0] = w_update0(0,is_obs[0], w, tau2, y,xb, neigh,U[0], B,  F_scale);
  
  for(int j=1 ; j < n ; j++){
    out[j] = w_update(j,is_obs[j], w, tau2, y,xb, neigh,U[j], B,  F_scale);
  } 
  
  return out;
}

//void bet_update(NumericVector &bet,vec mb, mat Vbinv,vec w,double tau2,vec y,mat X,mat XtX){ 
//  vec mu_new = Vbinv * mb + 1/tau2 * (X.t() * (y- w)) ;
//  mat V_inv = Vbinv + XtX/tau2 ;
//  mat V_new = inv(V_inv);
//  bet = my_mvrnorm(1,V_new*mu_new,V_new);
//}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]] 
vec bet_update(mat Vbinv,mat xtx, vec Vbinvmb,vec w,double tau2,vec y,mat X){ 
// int n = w.size();
//  int q = X.n_cols;
  vec mu_new = Vbinvmb + 1/tau2 * X.t() * (y-w) ;
  mat V_inv = Vbinv +xtx / tau2 ;
  mat V_new = V_inv.i() ;
  return my_mvrnorm(V_new * mu_new,V_new);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]] 

double tau2_update(double a, double b,vec w,vec y,vec xb){ 
  
  int n = w.size();
  double a_s = a + (double)n/2.0;
  double b_s = b + sum(pow(y - xb - w,2.0))/2.0;
  double prec = rgamma(1,a_s,1/b_s)[0];
  double tau2 = 1/prec;
  return tau2;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]] 

double sig2_update(double a, double b,List B,List F,vec w,List neigh){ 
  int n = w.size();
  NumericVector res(n);
  res[0] = w[0];
  // Rcpp::Rcout << "res = " << res[0] << std::endl;
  
  for(int i=1 ; i < n; i++){
    uvec t_neigh = neigh[i];
    vec B_temp = B[i];
    vec Bw =  B_temp.t() * w.elem(t_neigh) ;
    mat Ftemp = F[i];
    res[i] =  pow(w[i] - Bw[0],2.0) / Ftemp(0,0);
    // if( res[i] == NA ){
    //   Rcpp::Rcout << i << " is na | " << "  F is "<< Ftemp(0,0) << std::endl;
    // }
  }
  
  
  double prec = rgamma(1,n/2 + a,1/(b + sum(res)/2) )[0] ;
  double sig2 = 1 / prec ;
  return sig2;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]] 

List cov_update(vec covs, vec cand,List B,List F,double current_like,double sig2,
                int mod_choice,vec w,List neigh,int count, 
                List circ_day, List circ_week, List times){
  int n = w.size();
  
  /// covs, output &B,List &F &current_like, count
  
  List cors_temp = cor_list( circ_day,circ_week,times,cand,mod_choice,n );
  List Bs_temp = Bmats( cors_temp, n);
  List Fs_temp = Fmats( cors_temp, Bs_temp, n);
  double Ftemp = Fs_temp[0];
  double new_like = dnorm(w[0],0,sqrt(sig2 * Ftemp),true);
  //  Rcpp::Rcout << "new like = " << new_like << std::endl;
  for(int i=1 ; i < n; i++){
    uvec t_neigh = neigh[i];
    vec B_temp = Bs_temp[i];
    vec Bw =  B_temp.t() * w.elem(t_neigh) ;
    mat Ftemp = Fs_temp[i];
    new_like += dnorm(w[i],Bw[0],sqrt(sig2 * Ftemp(0,0)),true);
  }
  
  // Rcpp::Rcout << "old like = " << current_like << "new like = " << new_like << std::endl;
  
  double logu = log(runif(1,0,1)[0]) ;
  if(new_like - current_like > logu){
    covs = cand ;
    B = Bs_temp;
    F = Fs_temp;
    current_like = new_like;
    count = count + 1;
  }
  
  List res; res["covs"] = covs; res["B"] = B; res["F"] = F; 
  res["like"] = current_like; res["count"] = count;
  
  return res;
  
}

// List fit_NNGP(vec y, vec times,vec lon, vec lat,List neighs, 
//               int S, int burn,int mod_choice, int thin=1,int tune=100){
//   
//   int n = y.size();
//   int p = 4;
//   
//   int n_tot = S*thin + burn ;
//   double a_t = 1;
//   double b_t = 1;
//   double a_s = 1;
//   double b_s = 1;
//   //  double a_r = 0.1;
//   //  double b_r = 0.1;
//   //  vec mb(p,fill::zeros);
//   //  mat Vb= 1000.0 * eye<mat>(p,p);
//   //  double mb = 0;
//   //  double Vbinv = 1e-3;
//   NumericVector vars(p,0.001);
//   NumericVector count(p,0.0);
//   
//   NumericVector res_mu(S);
//   NumericVector res_sig2(S);
//   NumericVector res_tau2(S);
//   NumericMatrix res_w(S,n);
//   NumericMatrix res_cov(S,p);
//   
//   double sig2 = 1.0;
//   double tau2 = 1.0;
//   double cs = 0.5;
//   double ct = 0.5;
//   double alp = 1;
//   double bet = .5;
//   
//   double mu = mean(y); 
//   vec xb = mu*ones(n);
//   
//   //  vec bet(p,fill::zeros);
//   vec w(n) ;
//   for(int i =0; i < n ;i++){
//     w[i] = y[i] - mu + rnorm(1,0,5)[0];
//   }
//   
//   List U = inv_neigh(neighs, n);
//   List dists;
//   
//   if( (mod_choice == 3) | (mod_choice == 4) ){
//     dists = list_dists(lon,lat, neighs,2);
//   } else{
//     dists = list_dists(lon,lat, neighs,1);
//   } 
//   
//   List time_difs = list_time_difs(times,neighs);
//   List cors = cor_list( dists,time_difs,cs,ct,alp,bet,mod_choice,n );
//   List Bs = Bmats( cors, n);
//   List Fs = Fmats( cors, Bs, n); 
//   double Ftemp = Fs[0];
//   double current_like = dnorm(w[0],0,sqrt(sig2 * Ftemp),true);
//   
//   for(int i=1 ; i < n; i++){
//     uvec t_neigh = neighs[i];
//     vec B_temp = Bs[i];
//     vec Bw =  B_temp.t() * w.elem(t_neigh) ;
//     mat Ftemp = Fs[i];
//     current_like += dnorm(w[i],Bw[0],sqrt(sig2 * Ftemp(0,0)),true);
//   }
//   
//   for (int i = 1; i < n_tot; i++) {
//     
//     tau2_update(tau2,a_t , b_t, w,y,xb);
//     // Rcpp::Rcout << "tau2 = " << tau2 << std::endl;
//     sig2_update(sig2,a_s , b_s, Bs,Fs, w, neighs);
//     // Rcpp::Rcout << "sig2 = " << sig2 << std::endl;
//     cs_update(cs,ct,alp,bet,Bs,Fs,current_like, sig2, mod_choice, w, 
//               neighs, vars,count, dists, time_difs);
//     // Rcpp::Rcout << "cs = " << cs << std::endl;
//     ct_update(cs,ct,alp,bet,Bs,Fs,current_like, sig2, mod_choice, w, 
//               neighs, vars,count, dists, time_difs);
//     // Rcpp::Rcout << "ct = " << ct << std::endl;
//     alp_update(cs,ct,alp,bet,Bs,Fs,current_like, sig2, mod_choice, w, 
//                neighs, vars,count, dists, time_difs);
//     // Rcpp::Rcout << "alp = " << alp << std::endl;
//     bet_update(cs,ct,alp,bet,Bs,Fs,current_like, sig2, mod_choice, w, 
//                neighs, vars,count, dists, time_difs);
//     // Rcpp::Rcout << "bet = " << bet << std::endl;
//     
//     List F = Scale_F(Fs, sig2, n);
//     
//     w_update0(0, w, tau2, y,xb, neighs,U[0], Bs,  F);
//     for(int j=1 ; j < n ; j++){
//       w_update(j, w, tau2, y,xb, neighs,U[j], Bs,  F);
//     } 
//     double Ftemp = Fs[0];
//     current_like = 0;
//     current_like += dnorm(w[0],0,sqrt(sig2 * Ftemp),true);
//     
//     for(int j=1 ; j < n; j++){
//       uvec t_neigh = neighs[j];
//       vec B_temp = Bs[j];
//       vec Bw =  B_temp.t() * w.elem(t_neigh) ;
//       mat Ftemp = Fs[j];
//       current_like += dnorm(w[j],Bw[0],sqrt(sig2 * Ftemp(0,0)),true);
//     }
//     
//     if ( (i >= burn) & (i % thin == 0)) {
//       res_mu[(i - burn)/thin] = mu;
//       res_sig2[(i - burn)/thin] = sig2;
//       res_tau2[(i - burn)/thin] = tau2;
//       res_cov((i - burn)/thin,0) = cs;
//       res_cov((i - burn)/thin,1) = ct;
//       res_cov((i - burn)/thin,2) = alp;
//       res_cov((i - burn)/thin,3) = bet;
//       
//       for(int j=0; j < n ; j++){
//         res_w((i - burn)/thin,j) = w[j] ; 
//       }
//       
//     }
//     
//     
//     Rprintf("\r Iteration number: \v%d", i+1);
//     R_FlushConsole();
//     R_ProcessEvents();
//   }
//   
//   List res; res["mu"] = res_mu; res["sig2"] = res_sig2; res["w"] = res_w; res["tau2"] = res_tau2;
//   res["cov"] = res_cov ;res["vars"] = vars;
//   return(res);
// }

// List fit_NNGP(vec y, vec times,vec lon, vec lat,List neighs, 
//               int S, int burn,int mod_choice, int thin=1,int tune=100){

// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// NumericMatrix pred_NNGP(vec times,vec hours,vec grid_times,vec grid_hours,List grid_neighs, NumericVector mu, 
//                         NumericVector sig2, NumericVector tau2,NumericMatrix w,NumericMatrix cov_pars,int mod_choice){
  
//   int g = grid_times.size();
//   int rep = cov_pars.nrow();
  
//   NumericMatrix res(rep,g);
  
//   List circ_difs = list_grid_circ_difs(hours,grid_hours, grid_neighs);
//   List time_difs = list_grid_time_difs(times,grid_times,grid_neighs);
  
//   for(int i=0 ; i < rep; i++){
//     vec rho = cov_pars(i,_);
//     List cors = cor_list( circ_difs,time_difs,rho,mod_choice,g);
//     mat this_thing = cors[0];
//     List Bs = Bmats_g( cors, g);
//     List Fs = Fmats_g( cors, Bs, g);
//     vec w_sub = w(i,_);
    
//     for(int j=0 ; j < g; j++){
      
//       uvec t_neigh = grid_neighs[j];
//       vec B_temp = Bs[j];
//       vec temp = w_sub.elem(t_neigh);
//       vec Bw =  B_temp.t() * w_sub.elem(t_neigh) ;
//       mat Ftemp = Fs[j];
//       res(i,j) = mu[i] + rnorm(1,Bw[0],sqrt(sig2[i]*Ftemp(0,0) + tau2[i]))[0];
//     }
//     Rprintf("\r Prediction number: \v%d", i+1);
//     R_FlushConsole();
//     R_ProcessEvents();
//   }
//   return res;
// }
