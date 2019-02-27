library(foreign)
library(dplyr )
STATAfile <- read.dta("../data/complete.dta")

STATAfile$rac1<-as.numeric(STATAfile$rac1)-1
STATAfile$rac2<-as.numeric(STATAfile$rac2)-1
STATAfile$rac3<-as.numeric(STATAfile$rac3)-1

if (by_country==TRUE){
  STATAfile <- STATAfile[STATAfile$C==which_c,]
}

data.df  <-STATAfile[, c("imm1", "imm2", "imm3", "imm4",
                         "eco1", "eco2", "eco3", 
                         "com1", "com2", "com3", "com4",
                         "alt1", "alt2", "alt3",
                         "rac1", "rac2", "rac3", "rac4", "rac5",
                         "edu_high", "age", "unemployed", "retired",
                         "house_inc_decile", "minority", "male",
                         "dom_city", "dom_sub", "dom_town", "dom_vill", "dom_farm",
                         "country", "C")]
rm(STATAfile)
data.df <- data.df[complete.cases(data.df),]
if (reduced_data==TRUE){
  data.df <- sample_n(data.df, 200)
}
C  <- data.df$C
nC <- length(names(table(data.df$country)))


N <- dim(data.df)[1]

com <- matrix(0,4,N)
com[1,] <-  data.df$com1
com[2,] <-  data.df$com2
com[3,] <-  data.df$com3
com[4,] <-  data.df$com4

eco <- matrix(0,3,N)
eco[1,] <-  data.df$eco1
eco[2,] <-  data.df$eco2
eco[3,] <-  data.df$eco3

alt <- matrix(0,3,N)
alt[1,] <- round(data.df$alt1)
alt[2,] <- round(data.df$alt2)
alt[3,] <- round(data.df$alt3)

rac <- matrix(0,5,N)
rac[1,] <-  data.df$rac1
rac[2,] <-  data.df$rac2
rac[3,] <-  data.df$rac3
rac[4,] <-  data.df$rac4
rac[5,] <-  data.df$rac5

imm <- matrix(0,4,N)
imm[1,] <-  data.df$imm1
imm[2,] <-  data.df$imm2
imm[3,] <-  data.df$imm3
imm[4,] <-  data.df$imm4

cov <-  matrix(0,N, 11)
cov[,1]  <- data.df$edu_high -mean(data.df$edu_high)
cov[,2]  <- (data.df$age     -mean(data.df$age))/(2*sd(data.df$age))
cov[,3]  <- (as.numeric(data.df$house_inc_decile)-mean(as.numeric(data.df$house_inc_decile)))/(2*sd(as.numeric(data.df$house_inc_decile)))
cov[,4]  <- data.df$unemployed-mean(data.df$unemployed)
cov[,5]  <- data.df$retired  -mean(data.df$retired)
cov[,6]  <- data.df$male     -mean(data.df$male)
cov[,7]  <- data.df$minority -mean(data.df$minority)
cov[,8]  <- data.df$dom_sub  -mean(data.df$dom_sub  )
cov[,9]  <- data.df$dom_town -mean(data.df$dom_town )
cov[,10] <- data.df$dom_vill -mean(data.df$dom_vill )
cov[,11] <- data.df$dom_farm -mean(data.df$dom_farm )

weight <-   data.df$pweight

lkj_const <- 10

K <- dim(cov)[2]

