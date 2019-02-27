rm(list=ls())
library(rstan)
library(foreign)
library(dplyr )
library(magrittr)

data.df <- 
  read.dta("../raw_data/ESS7e02_2.dta") %>% 
  select(cntry,
         imsmetn,
         imdfetn,
         eimpcnt,
         impcntr,
         imtcjob,
         imbleco,
         imbgeco,
         imueclt,
         rlgueim,
         imwbcrm,
         imwbcnt,
         iphlppl,
         ipeqopt,
         ipudrst,
         hlpfmly,
         smegbli,
         smegbhw,
         smctmbe,
         qfimwht,
         imdetbs,
         imdetmr,
         agea,
         gndr,
         edulvlb,
         blgetmg,
         hinctnta,
         hincsrca,
         uempla,
         domicil) %>% 
  filter(cntry!="AT",
         cntry!="EE",
         cntry!="CZ") %>% 
  droplevels() %>% 
  transmute(
    country = cntry,
    C = as.numeric(as.factor(country)),
    imm1 = (4-as.numeric(imsmetn))+1,
    imm2 = (4-as.numeric(imdfetn))+1,
    imm3 = (4-as.numeric(eimpcnt))+1,
    imm4 = (4-as.numeric(impcntr))+1,
    eco1 = as.numeric(imtcjob)+1,
    eco2 = as.numeric(imbleco)+1,
    eco3 = as.numeric(imbgeco)+1,
    com1 = as.numeric(imueclt)+1,
    com2 = as.numeric(rlgueim)+1,
    com3 = as.numeric(imwbcrm)+1,
    com4 = as.numeric(imwbcnt)+1,
    alt1 = 7-as.numeric(iphlppl),
    alt2 = 7-as.numeric(ipeqopt),
    alt3 = 7-as.numeric(ipudrst),
    alt4 = (-as.numeric(hlpfmly)+2)+1,
    rac1 = 3-as.numeric(smegbli),
    rac2 = 3-as.numeric(smegbhw),
    rac3 = 3-as.numeric(smctmbe),
    rac4 = as.numeric(qfimwht)+1,
    rac5 = as.numeric(imdetbs)+1,
    rac6 = as.numeric(imdetmr)+1,
    age = agea,
    male = as.numeric(gndr=="Male"),
    edu_high = as.numeric(as.numeric(edulvlb)>=21 & edulvlb!="Other"),
    minority = as.numeric(blgetmg == "Yes"),
    house_inc_decile = as.numeric(hinctnta),
    retired = as.numeric(hincsrca == "Pensions"),
    unemployed = as.numeric(uempla == "Marked"),
    dom_city = as.numeric(domicil == "A big city"),
    dom_sub = as.numeric(domicil == "Suburbs or outskirts of big city"),
    dom_town = as.numeric(domicil == "Town or small city"),
    dom_vill = as.numeric(domicil == "Country village"),
    dom_farm = as.numeric(domicil == "Farm or home in countryside")
    ) %>%
  filter(complete.cases(.)) 
  
save_items <- 
  function(dataframe, filename){
   C  <- dataframe$C
   nC <- length(names(table(dataframe$country)))
   N <- nrow(dataframe)
   
   com <- matrix(0,4,N)
   com[1,] <-  dataframe$com1
   com[2,] <-  dataframe$com2
   com[3,] <-  dataframe$com3
   com[4,] <-  dataframe$com4
   
   eco <- matrix(0,3,N)
   eco[1,] <-  dataframe$eco1
   eco[2,] <-  dataframe$eco2
   eco[3,] <-  dataframe$eco3
   
   alt <- matrix(0,3,N)
   alt[1,] <- round(dataframe$alt1)
   alt[2,] <- round(dataframe$alt2)
   alt[3,] <- round(dataframe$alt3)
   
   rac <- matrix(0,5,N)
   rac[1,] <-  dataframe$rac1
   rac[2,] <-  dataframe$rac2
   rac[3,] <-  dataframe$rac3
   rac[4,] <-  dataframe$rac4
   rac[5,] <-  dataframe$rac5
   
   imm <- matrix(0,4,N)
   imm[1,] <-  dataframe$imm1
   imm[2,] <-  dataframe$imm2
   imm[3,] <-  dataframe$imm3
   imm[4,] <-  dataframe$imm4
   
   cov <-  matrix(0,N, 11)
   cov[,1]  <- dataframe$edu_high -mean(dataframe$edu_high)
   cov[,2]  <- (dataframe$age     -mean(dataframe$age))/(2*sd(dataframe$age))
   cov[,3]  <- {(as.numeric(dataframe$house_inc_decile)
                 -mean(as.numeric(dataframe$house_inc_decile)))/
       (2*sd(as.numeric(dataframe$house_inc_decile)))
   }
   
   cov[,4]  <- dataframe$unemployed-mean(dataframe$unemployed)
   cov[,5]  <- dataframe$retired  -mean(dataframe$retired)
   cov[,6]  <- dataframe$male     -mean(dataframe$male)
   cov[,7]  <- dataframe$minority -mean(dataframe$minority)
   cov[,8]  <- dataframe$dom_sub  -mean(dataframe$dom_sub  )
   cov[,9]  <- dataframe$dom_town -mean(dataframe$dom_town )
   cov[,10] <- dataframe$dom_vill -mean(dataframe$dom_vill )
   cov[,11] <- dataframe$dom_farm -mean(dataframe$dom_farm )
   
   K <- dim(cov)[2]
   
   lkj_const <- 10
   
   stan_rdump(
     c(
       'N',
       'K',
       'C',
       'nC',
       'com',
       'eco',
       'alt',
       'rac',
       'imm',
       'cov',
       'lkj_const'), 
     file=paste0("../processed_data/",filename,".Rdump"))
  }

save_items(data.df, "dat")
