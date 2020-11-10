library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(xtable)
library(Matrix)
library(GpGp)
library(inline)
library(MASS)
#RcppArmadillo.package.skeleton()

try(setwd("C:/Users/philaw/Box/Research/Porcu/Torus/code/Space_circle/cloud_cover_w"),silent=TRUE)
rm(list=ls())

########################## Load Data

#source("temp_dat.R")
load("dat_prep_full.RData")
source("fit_nngp.R")

Rcpp::sourceCpp("nngp_circle.cpp")
y = dat_ordered$tcdc

in_function = function(x,val,alp=0.1){
  CI = quantile(x,c(alp/2,1 - alp/2))
  1*(val <= CI[2] & val >= CI[1])
}

mod_choice = 8

NN = NNarray[,-1]

S = 20e3
burn = 10e3

mod_out = fit_NNGP(dat_ordered$tcdc, NN, dat_ordered$time,
                   dat_ordered$lon, dat_ordered$lat, 
                   S = S, is_obs = rep(TRUE,n), burn = burn, 
                   mod_choice = mod_choice, thin=1, tune=100)

save.image("mod8.RData")
