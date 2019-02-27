library(rstan)
library(BH)
library(ggplot2)
library(shinystan)
library(plyr)
library(rstan)
library(reshape2)
library(xtable)
library(foreign)


reduced_data = FALSE
by_country   = FALSE
which_c      = 1

source("dat.r"                       , echo=TRUE)
stan_rdump(
  c('N','K','C','nC','com','eco','alt', 'rac','imm','cov', 'lkj_const'), 
  file="dat.Rdump")


for (i in 1:18){
  by_country   <- TRUE
  which_c      <- i
  source("dat.r"                       , echo=TRUE)
  stan_rdump(
    c('N','K','C','nC','com','eco','alt', 'rac','imm','cov','lkj_const'),
    file = paste("dat_",which_c,".Rdump", sep=""))
}



