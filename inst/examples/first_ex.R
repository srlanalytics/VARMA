library(devtools)
library(data.table)
library(Rcpp)
library(parallel)
load_all("~/bdfm")
load_all("~/seasonality")
source("~/OttoQuant_Inc/R_backend/R/general_utils.R")
load("C:/Users/seton/Dropbox/data/gKNi_raw.RData")
source("~/OttoQuant_Inc/R_backend/R/process_inputs.R")

targets <- c("bloomberg.imports.us", "bloomberg.exports.us")
as_of <- as.Date("2018-01-01")
y_is_SA <- T
y <- call_data_wide(targets, dt = dt_raw, as_of = as_of, series_id = "series_name", pub_date = "pub_time")

processed <- process_targets(y, y_is_SA)

Y <- processed$data[,-1,with = F]

ts.plot(Y)

library(MTS)

Y <- Y[-1,]

pred <- VARMA(Y, p = 3, q = 3, include.mean = F)

phi=matrix(c(0.2,-0.6,0.3,1.1),2,2); theta=matrix(c(-0.5,0,0,-0.5),2,2)
sigma=diag(2)
m1=VARMAsim(300,arlags=c(1),malags=c(1),phi=phi,theta=theta,sigma=sigma)
zt=m1$series
m2=VARMA(zt,p=1,q=1,include.mean=FALSE)