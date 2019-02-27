library(tidyr)
library(xtable)
library(dplyr)
library(tibble)
library(magrittr)
library(rstan)
options(mc.cores = parallel::detectCores())
library(shinystan)
library(foreign)
library(countrycode)
library(doParallel)
library(stargazer)
rm(list=ls())

starTeX <- function(filename, repl, ...) {
  tmp.prnt <- capture.output(stargazer(..., type = "latex"))
  for (i in 1:nrow(repl)) {
    tmp.prnt <- gsub(repl[i, 1], repl[i, 2], tmp.prnt, fixed = T)
  }
  write(tmp.prnt, file = filename)
}


STATAfile <- read.dta("../data/complete.dta")
data.df  <-STATAfile[, c("edu_high","country", "C")]
rm(STATAfile)
data.df <- data.df[complete.cases(data.df),]
countrycode <- names(table(data.df$country))
countryname <- 
  countrycode(countrycode, origin="iso2c", destination = "country.name")

folder <- "M:/mcmc_stash/cmdstan18.0/output_by_c_"

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

countrylist <-  c(1:18)

cl <- makeCluster(6)
registerDoParallel(cl)
post <-
  foreach(j = countrylist,
          .packages = c("rstan", "shinystan")
          ) %dopar% {
    as.shinystan(
      read_stan_csv(
        c(
          paste0(folder, j, "_chain_1.csv"),
          paste0(folder, j, "_chain_2.csv"),
          paste0(folder, j, "_chain_3.csv"),
          paste0(folder, j, "_chain_4.csv")
        )
      )
    )
}
stopCluster(cl)

diagnostics <- data.frame(
  country   = countrylist,
  rhat      = rep(NA, length(countrylist)),
  diverging = rep(NA, length(countrylist))
  )

for(j in 1:length(countrylist)){
  diagnostics$diverging[j] <-
    sum(retrieve(post[[j]], "diverging"))

  diagnostics$rhat[j] <-
    retrieve(post[[j]], "rhats") %>%
    as.data.frame() %>%
    dplyr::arrange(., -.) %>%
    head(n = 1)
}

pick.estimates <- function(X,name){
  X[grep(name, names(X), fixed=T)]
}

diagnostics

post.means <- lapply(post, retrieve, "means")
post.sd    <- lapply(post, retrieve, "sd")

for(i in 1:18){
  edu.tmp <-
    data.df %>%
    filter(country==names(table(country))[i]) %>% 
    select(edu_high)
  
  dt      <- density(pick.estimates(post.means[[i]], "imm_lat"))
  d.edu   <- density(pick.estimates(post.means[[i]], "imm_lat")[edu.tmp==TRUE ])
  d.noedu <- density(pick.estimates(post.means[[i]], "imm_lat")[edu.tmp==FALSE])
  
  yrange  <- range(dt$y, d.noedu$y, d.edu$y)
  
  pdf(paste0("output/imm_pref_",countrynames[i],".pdf"))
    plot(dt,xlab = "Latent immigration preferences",
      main = countrynames[i],
      ylim=yrange)
    lines(d.edu, col="blue")
    lines(d.noedu, col="red")
    legend("topleft",
           c("Total", "High Education", "Low Education"),
           col=c("black", "blue", "red"),
           lty=1,
           bty='n')
  dev.off()
}

add.paranteses <- function(x){
  paste0("(",x,")")
  }

coef.means <- 
  sapply(post.means, pick.estimates, "beta") %>%
  rbind(sapply(post.means, pick.estimates, "delta[")) %>%
  round(digits=2)
coef.sd    <- 
  sapply(post.sd   , pick.estimates, "beta")%>%
  rbind(sapply(post.sd, pick.estimates, "delta[")) %>%
  round(digits=2) %>%
  apply(2, add.paranteses)
rownames(coef.sd) <- paste0("sd_" ,rownames(coef.means))
coef <- 
  rbind(coef.means, coef.sd) %>%
  .[order(rep((1:(nrow(.)/2))*10,2)+rep(1:2, each=(nrow(.)/2))),]

colnames(coef) <- countrycode


repl <- data.frame(
  what =  
    c(
      paste0("sd\\_beta[", 1:4, "]"),
      paste0("sd\\_delta[", 1:11, "]"),
      paste0("beta[", 1:4, "]"),
      paste0("delta[", 1:11, "]"),
      "Edu",
      "cccccccccc",
      "\\begin{tabular}","\\end{tabular}",
      "\\hline","\\\\[-1.8ex]",
      "ccc}",
      "HU \\\\" ,
      "SI \\\\" ,
      "\\end{longtable}",
      "\\bottomrule"
    ),
  with = 
    c(
      rep("",4),
      rep("",11),
      paste0("$\\beta_{", c("Eco", "Com", "Rac", "Alt"), "}$"),
      deltanames,
      "\\midrule \\ Edu",
      "lccccccccc",
      "\\begin{longtable}","\\end{longtable}",
      "","",
      "ccc} \\toprule",
      "HU \\\\ \\midrule",
      "SI \\\\ \\midrule",
      "\\bottomrule \\end{longtable}",
      "\\bottomrule
      \\caption*{\\parbox[t]{17cm}{\\scriptsize{{\\it Note: \\textnormal{
      The table shows the estimates of the model explained in section 
      \\ref{sec:model}, where we use four survey questions (IMM1, IMM2, IMM3 
      and IMM4) to identify individual, latent preference for immigration. 
      Plain numbers indicate posterior means and numbers in parentheses show 
      thestandard deviation of the posterior distribution. The model is 
      estimated on each country separately. Note that the variables age and 
      income decile have been standardized by two times their standard 
      deviation. 
      }}}}}"
    )
    )


starTeX("output/beta1.tex", repl, coef[,1:9])
starTeX("output/beta2.tex", repl, coef[,10:18])

