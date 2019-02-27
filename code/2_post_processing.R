library(tidyr)
library(xtable)
library(dplyr)
library(tibble)
library(magrittr)
library(rstan)
options(mc.cores = parallel::detectCores())
library(shinystan)
rm(list=ls())


folder <- "M:/mcmc_stash/cmdstan17.0 - copy/output_"
#folder <- "M:/mcmc_stash/cmdstan/output_"
#folder <- "cmdstan16.0/output_"
#folder <- "M:/mcmc_stash/cmdstan_small/output_"
#folder <- "C:/cmdstan/cmdstan-2.17.1/cmdstan17.0/output_"

deltanames <- c(
  "Edu",
  "Age",
  "Inc.Decile",
  "Unemployed",
  "Retired",
  "Male",
  "Minority",
  "Suburb",
  "Town",
  "Village",
  "Farm"
  )

rhonames <- c( 
  "rho_eco_com", 
  "rho_eco_rac", 
  "rho_eco_alt", 
  "rho_com_rac", 
  "rho_com_alt", 
  "rho_rac_alt"
  )

post_multi<-
  as.shinystan(
    read_stan_csv(
      c(
        paste0(folder, "multi", "_1.csv"),
        paste0(folder, "multi", "_2.csv"),
        paste0(folder, "multi", "_3.csv"),
        paste0(folder, "multi", "_4.csv")
      )
    )
  )

save(post_multi, "post.Rdata")
# Diverging: 
sum(retrieve(post_multi,"diverging"))

# Rhat > 1.05: 
retrieve(post_multi,"rhats") %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  filter(.>1.05) %>%
  dplyr::arrange(.,-.) %>%
  head()

launch_shinystan(post_multi)

s_multi  <- post_multi @summary %>%
  head(.,-1) %>%
  as.data.frame %>%
  select("mean", "sd") %>%
  rownames_to_column() %>%
  slice(-grep("lat_f", .$rowname))

lat_means  <- 
  post_multi @summary %>%
  head(.,-1) %>%
  as.data.frame %>%
  select("mean") %>%
  rownames_to_column() %>%
  slice(
    c(grep("imm_lat", .$rowname),
      grep("lat_f"  , .$rowname)
      )) %>%
  slice(-grep("sd", .$rowname)) %>%
  mutate(
    n = c(1:(nrow(.)/5),
          rep(1:(nrow(.)/5), each=4)
          ),
    v = c(rep(1, nrow(.)/5),
          rep(2:5,nrow(.)/5)
          )) %>%
  select(mean,n,v) %>%
  spread(key = c("v"), value="mean")
  
colnames(lat_means) <- c("n","imm", "eco","com","rac","alt")


complete.multi   <-as.data.frame(matrix(NA,2*(length(deltanames)+4), 5))

colnames(complete.multi) <- c("Imm","Eco","Com","Rac","Alt")
i.delta.e <- grep("delta_e", s_multi$rowname )
i.delta.a <- grep("delta_a", s_multi$rowname )
i.delta.c <- grep("delta_c", s_multi$rowname )
i.delta.r <- grep("delta_r", s_multi$rowname )
i.delta   <- grep("delta[", s_multi$rowname, fixed = T)

s_multi[,2]<- round(s_multi[,2],digits=3)
s_multi[,3]<- round(s_multi[,3],digits=3)

complete.multi[seq(from=1, to=2*length(deltanames), by=2)+8,]<-
 cbind(
  s_multi[i.delta  ,]$mean,
  s_multi[i.delta.e,]$mean,
  s_multi[i.delta.c,]$mean,
  s_multi[i.delta.r,]$mean,
  s_multi[i.delta.a,]$mean
  )

complete.multi[seq(from=2, to=2*length(deltanames), by=2)+8,]<-
  cbind(
    s_multi[i.delta  ,]$sd,
    s_multi[i.delta.e,]$sd,
    s_multi[i.delta.c,]$sd,
    s_multi[i.delta.r,]$sd,
    s_multi[i.delta.a,]$sd
  )

complete.multi[seq(from=1,to=8,by=2),1]<-
  s_multi[grep("beta",s_multi$rowname),]$mean

complete.multi[seq(from=2,to=8,by=2),1]<-
  s_multi[grep("beta",s_multi$rowname),]$sd

complete.multi$Variable     <- rep("  ", dim(complete.multi)[1])
complete.multi$Variable[seq(1, dim(complete.multi)[1], 2)]  <- 
  c("Eco","Com","Rac","Alt",deltanames)
complete.multi<-complete.multi[,c(6,1:5)]

print(xtable(complete.multi),include.rownames = FALSE)

smpls <- grep("beta",names(post_multi@posterior_sample[1,1,])) %>%
  post_multi@posterior_sample[,,.] %>%
  plyr::adply(., 2)

plot(density(lat_means$imm),
     main = "Latent preferences for immigration")
