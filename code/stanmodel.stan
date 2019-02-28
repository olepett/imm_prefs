functions{

// Function for calculating median of a ordered vector. If even number of 
// entries, return mean of two entries in center. 
real median(vector ordered_vector){
  int midpoint;
  midpoint = rows(ordered_vector)/2;
  if(2*midpoint<rows(ordered_vector)){
    return ordered_vector[midpoint+1];
  }
  else{
    return (ordered_vector[midpoint]+ordered_vector[midpoint+1])/2;
  }
}

// Function for calculting the cumulative sum of a vector, and then 
// subtracting the median from the result. 
vector de_median(vector par_vector){
  return cumulative_sum(par_vector)-median(cumulative_sum(par_vector));
}

}


data{
real<lower=0> lkj_const; // Shrinkage on LKJ-prior
int<lower=0> N;         // no. obs. 
int<lower=0> K;         // no. covariates
int C[N];               // vector with indices of country from each obs. 
int nC;                 // Total no. of countries
int com[4,N];           // Responses along com. dimension   
int eco[3,N];           // Responses along eco. dimension   
int alt[3,N];           // Responses along alt. dimension   
int rac[5,N];           // Responses along rac. dimension   
int imm[4,N];           // Responses along imm. dimension   
matrix[N,K] cov;        // Covariates   
}

parameters{
simplex[10] theta_com1;      // Relative distance of cutpoints
simplex[10] theta_com2;      // Relative distance of cutpoints
simplex[10] theta_com3;      // Relative distance of cutpoints
simplex[10] theta_com4;      // Relative distance of cutpoints
simplex[10] theta_eco1;      // Relative distance of cutpoints
simplex[10] theta_eco2;      // Relative distance of cutpoints
simplex[10] theta_eco3;      // Relative distance of cutpoints
simplex[5]  theta_alt1;      // Relative distance of cutpoints
simplex[5]  theta_alt2;      // Relative distance of cutpoints
simplex[5]  theta_alt3;      // Relative distance of cutpoints
simplex[10] theta_rac4;      // Relative distance of cutpoints
simplex[10] theta_rac5;      // Relative distance of cutpoints
simplex[3]  theta_imm1;      // Relative distance of cutpoints
simplex[3]  theta_imm2;      // Relative distance of cutpoints
simplex[3]  theta_imm3;      // Relative distance of cutpoints
simplex[3]  theta_imm4;      // Relative distance of cutpoints
real<lower=0> kappa_com1;    // Spread of cutpoints
real<lower=0> kappa_com2;    // Spread of cutpoints
real<lower=0> kappa_com3;    // Spread of cutpoints
real<lower=0> kappa_com4;    // Spread of cutpoints
real<lower=0> kappa_eco1;    // Spread of cutpoints
real<lower=0> kappa_eco2;    // Spread of cutpoints
real<lower=0> kappa_eco3;    // Spread of cutpoints
real<lower=0> kappa_alt1;    // Spread of cutpoints
real<lower=0> kappa_alt2;    // Spread of cutpoints
real<lower=0> kappa_alt3;    // Spread of cutpoints
real<lower=0> kappa_rac4;    // Spread of cutpoints
real<lower=0> kappa_rac5;    // Spread of cutpoints
real<lower=0> kappa_imm1;    // Spread of cutpoints
real<lower=0> kappa_imm2;    // Spread of cutpoints
real<lower=0> kappa_imm3;    // Spread of cutpoints
real<lower=0> kappa_imm4;    // Spread of cutpoints
real cutoffs_rac1;           // Overall mean of cutpoints
real cutoffs_rac2;           // Overall mean of cutpoints
real cutoffs_rac3;           // Overall mean of cutpoints
real mu_com1;                // Overall mean of cutpoints
real mu_com2;                // Overall mean of cutpoints
real mu_com3;                // Overall mean of cutpoints
real mu_com4;                // Overall mean of cutpoints
real mu_eco1;                // Overall mean of cutpoints
real mu_eco2;                // Overall mean of cutpoints
real mu_eco3;                // Overall mean of cutpoints
real mu_alt1;                // Overall mean of cutpoints
real mu_alt2;                // Overall mean of cutpoints
real mu_alt3;                // Overall mean of cutpoints
real mu_rac4;                // Overall mean of cutpoints
real mu_rac5;                // Overall mean of cutpoints
real mu_imm1;                // Overall mean of cutpoints
real mu_imm2;                // Overall mean of cutpoints
real mu_imm3;                // Overall mean of cutpoints
real mu_imm4;                // Overall mean of cutpoints


real<lower=0> D_rac1;                // Discrimination of item
real<lower=0> D_rac2;                // Discrimination of item
real<lower=0> D_rac3;                // Discrimination of item
real<lower=0> D_com1;                // Discrimination of item
real<lower=0> D_com2;                // Discrimination of item
real<lower=0> D_com3;                // Discrimination of item
real<lower=0> D_com4;                // Discrimination of item
real<lower=0> D_eco1;                // Discrimination of item
real<lower=0> D_eco2;                // Discrimination of item
real<lower=0> D_eco3;                // Discrimination of item
real<lower=0> D_alt1;                // Discrimination of item
real<lower=0> D_alt2;                // Discrimination of item
real<lower=0> D_alt3;                // Discrimination of item
real<lower=0> D_rac4;                // Discrimination of item
real<lower=0> D_rac5;                // Discrimination of item
real<lower=0> D_imm1;                // Discrimination of item
real<lower=0> D_imm2;                // Discrimination of item
real<lower=0> D_imm3;                // Discrimination of item
real<lower=0> D_imm4;                // Discrimination of item


vector[4] beta;              // Coeffs. on latent factors
vector[K] delta;             // coeffs. on covariates, imm dimension
vector[K] delta_e;           // coeffs. on covariates, eco dimension
vector[K] delta_c;           // coeffs. on covariates, com dimension
vector[K] delta_r;           // coeffs. on covariates, rac dimension
vector[K] delta_a;           // coeffs. on covariates, alt dimension
vector[nC-1] cfe_raw;        // country-fixed effects, imm dimension
vector[nC-1] cfe_raw_e;      // country-fixed effects, eco dimension
vector[nC-1] cfe_raw_c;      // country-fixed effects, com dimension
vector[nC-1] cfe_raw_r;      // country-fixed effects, rac dimension
vector[nC-1] cfe_raw_a;      // country-fixed effects, alt dimension

real<lower=0> spread_var;    // Hyper-par. for spread of cutpoints
real<lower=0> spread_mean;   // Hyper-par. for mean-cutpoints

vector[N] imm_lat;           // Latent pref. for immigration        
real<lower=0> sd_imm;        // Scale of imm_lat

vector[4] lat_f[N];          // Latent factors for eco, com, rac and alt
cholesky_factor_corr[4] L;   // Chol. fact. corr. for lat_f
vector<lower=0>[4] sd_lat_f; // Scale vec. for lat.f 

}


transformed parameters{
vector[nC] cfe;     // country-fixed effects, imm dimension
vector[nC] cfe_e;   // country-fixed effects, eco dimension
vector[nC] cfe_c;   // country-fixed effects, com dimension
vector[nC] cfe_r;   // country-fixed effects, rac dimension
vector[nC] cfe_a;   // country-fixed effects, alt dimension


ordered[10] cutoffs_com1; // Cutpoints
ordered[10] cutoffs_com2; // Cutpoints
ordered[10] cutoffs_com3; // Cutpoints
ordered[10] cutoffs_com4; // Cutpoints
ordered[10] cutoffs_eco1; // Cutpoints
ordered[10] cutoffs_eco2; // Cutpoints
ordered[10] cutoffs_eco3; // Cutpoints
ordered[5]  cutoffs_alt1; // Cutpoints
ordered[5]  cutoffs_alt2; // Cutpoints
ordered[5]  cutoffs_alt3; // Cutpoints
ordered[10] cutoffs_rac4; // Cutpoints
ordered[10] cutoffs_rac5; // Cutpoints
ordered[3]  cutoffs_imm1; // Cutpoints
ordered[3]  cutoffs_imm2; // Cutpoints
ordered[3]  cutoffs_imm3; // Cutpoints
ordered[3]  cutoffs_imm4; // Cutpoints

cutoffs_com1=de_median(theta_com1)*kappa_com1+mu_com1;
cutoffs_com2=de_median(theta_com2)*kappa_com2+mu_com2;
cutoffs_com3=de_median(theta_com3)*kappa_com3+mu_com3;
cutoffs_com4=de_median(theta_com4)*kappa_com4+mu_com4;
cutoffs_eco1=de_median(theta_eco1)*kappa_eco1+mu_eco1;
cutoffs_eco2=de_median(theta_eco2)*kappa_eco2+mu_eco2;
cutoffs_eco3=de_median(theta_eco3)*kappa_eco3+mu_eco3;
cutoffs_alt1=de_median(theta_alt1)*kappa_alt1+mu_alt1;
cutoffs_alt2=de_median(theta_alt2)*kappa_alt2+mu_alt2;
cutoffs_alt3=de_median(theta_alt3)*kappa_alt3+mu_alt3;
cutoffs_rac4=de_median(theta_rac4)*kappa_rac4+mu_rac4;
cutoffs_rac5=de_median(theta_rac5)*kappa_rac5+mu_rac5;
cutoffs_imm1=de_median(theta_imm1)*kappa_imm1+mu_imm1;
cutoffs_imm2=de_median(theta_imm2)*kappa_imm2+mu_imm2;
cutoffs_imm3=de_median(theta_imm3)*kappa_imm3+mu_imm3;
cutoffs_imm4=de_median(theta_imm4)*kappa_imm4+mu_imm4;

cfe[1] = 0;
cfe_e[1] = 0;
cfe_c[1] = 0;
cfe_r[1] = 0;
cfe_a[1] = 0;
for (i in 2:nC){
cfe[i] = cfe_raw[i-1];
cfe_e[i] = cfe_raw_e[i-1];
cfe_c[i] = cfe_raw_c[i-1];
cfe_r[i] = cfe_raw_r[i-1];
cfe_a[i] = cfe_raw_a[i-1];
}

}

model{
matrix[4, 4] L_Sigma;
matrix[N,4] lat_f_m;
D_rac1  ~lognormal(.5,1);
D_rac2  ~lognormal(.5,1);
D_rac3  ~lognormal(.5,1);
D_com1  ~lognormal(.5,1);
D_com2  ~lognormal(.5,1);
D_com3  ~lognormal(.5,1);
D_com4  ~lognormal(.5,1);
D_eco1  ~lognormal(.5,1);
D_eco2  ~lognormal(.5,1);
D_eco3  ~lognormal(.5,1);
D_alt1  ~lognormal(.5,1);
D_alt2  ~lognormal(.5,1);
D_alt3  ~lognormal(.5,1);
D_rac4  ~lognormal(.5,1);
D_rac5  ~lognormal(.5,1);
D_imm1  ~lognormal(.5,1);
D_imm2  ~lognormal(.5,1);
D_imm3  ~lognormal(.5,1);
D_imm4  ~lognormal(.5,1);
theta_com1~dirichlet(rep_vector(2,10));
theta_com2~dirichlet(rep_vector(2,10));
theta_com3~dirichlet(rep_vector(2,10));
theta_com4~dirichlet(rep_vector(2,10));
theta_eco1~dirichlet(rep_vector(2,10));
theta_eco2~dirichlet(rep_vector(2,10));
theta_eco3~dirichlet(rep_vector(2,10));
theta_alt1~dirichlet(rep_vector(2, 5));
theta_alt2~dirichlet(rep_vector(2, 5));
theta_alt3~dirichlet(rep_vector(2, 5));
theta_rac4~dirichlet(rep_vector(2,10));
theta_rac5~dirichlet(rep_vector(2,10));
theta_imm1~dirichlet(rep_vector(2, 3));
theta_imm2~dirichlet(rep_vector(2, 3));
theta_imm3~dirichlet(rep_vector(2, 3));
theta_imm4~dirichlet(rep_vector(2, 3));

spread_var~gamma(2,1);
kappa_com1~gamma(2, 1.0/spread_var);
kappa_com2~gamma(2, 1.0/spread_var);
kappa_com3~gamma(2, 1.0/spread_var);
kappa_com4~gamma(2, 1.0/spread_var);
kappa_eco1~gamma(2, 1.0/spread_var);
kappa_eco2~gamma(2, 1.0/spread_var);
kappa_eco3~gamma(2, 1.0/spread_var);
kappa_alt1~gamma(2, 1.0/spread_var);
kappa_alt2~gamma(2, 1.0/spread_var);
kappa_alt3~gamma(2, 1.0/spread_var);
kappa_rac4~gamma(2, 1.0/spread_var);
kappa_rac5~gamma(2, 1.0/spread_var);
kappa_imm1~gamma(2, 1.0/spread_var);
kappa_imm4~gamma(2, 1.0/spread_var);
kappa_imm3~gamma(2, 1.0/spread_var);
kappa_imm4~gamma(2, 1.0/spread_var);

spread_mean~gamma(2,1.0/2);
mu_com1~normal(0,spread_mean);
mu_com2~normal(0,spread_mean);
mu_com3~normal(0,spread_mean);
mu_com4~normal(0,spread_mean);
mu_eco1~normal(0,spread_mean);
mu_eco2~normal(0,spread_mean);
mu_eco3~normal(0,spread_mean);
mu_alt1~normal(0,spread_mean);
mu_alt2~normal(0,spread_mean);
mu_alt3~normal(0,spread_mean);
mu_rac4~normal(0,spread_mean);
mu_rac5~normal(0,spread_mean);
mu_imm1~normal(0,spread_mean);
mu_imm2~normal(0,spread_mean);
mu_imm3~normal(0,spread_mean);
mu_imm4~normal(0,spread_mean);
cutoffs_rac1~normal(0,spread_mean);
cutoffs_rac2~normal(0,spread_mean);
cutoffs_rac3~normal(0,spread_mean);

beta~normal(0,1);

delta~normal(0,1);
delta_e~normal(0,1);
delta_c~normal(0,1);
delta_r~normal(0,1);
delta_a~normal(0,1);

cfe_raw ~normal(0,1);
cfe_raw_a ~normal(0,1);
cfe_raw_c ~normal(0,1);
cfe_raw_r ~normal(0,1);
cfe_raw_a ~normal(0,1);

L~lkj_corr_cholesky(lkj_const);
sd_lat_f ~gamma(2, 1);
L_Sigma = diag_pre_multiply(sd_lat_f, L);

for (n in 1:N){
vector[4] muvec;
muvec[1] = cov[n,]*delta_e+cfe_e[C[n]];
muvec[2] = cov[n,]*delta_c+cfe_c[C[n]];
muvec[3] = cov[n,]*delta_r+cfe_r[C[n]];
muvec[4] = cov[n,]*delta_a+cfe_a[C[n]];
lat_f[n] ~ multi_normal_cholesky(muvec,L_Sigma);
}

for (n in 1:N){
lat_f_m[n,] = to_row_vector(lat_f[n]);
}

imm_lat   ~normal(
+beta[1]*lat_f_m[,1]/(2*sd(lat_f[,1]))
+beta[2]*lat_f_m[,2]/(2*sd(lat_f[,2]))
+beta[3]*lat_f_m[,3]/(2*sd(lat_f[,3]))
+beta[4]*lat_f_m[,4]/(2*sd(lat_f[,4]))
+cov*delta+cfe[C],sd_imm);

for (n in 1:N){
com[1,n]~ordered_logistic(D_com1 * lat_f[n,2] ,D_com1 * cutoffs_com1 );
com[2,n]~ordered_logistic(D_com2 * lat_f[n,2] ,D_com2 * cutoffs_com2 );
com[3,n]~ordered_logistic(D_com3 * lat_f[n,2] ,D_com3 * cutoffs_com3 );
com[4,n]~ordered_logistic(D_com4 * lat_f[n,2] ,D_com4 * cutoffs_com4 );
eco[1,n]~ordered_logistic(D_eco1 * lat_f[n,1] ,D_eco1 * cutoffs_eco1 );
eco[2,n]~ordered_logistic(D_eco2 * lat_f[n,1] ,D_eco2 * cutoffs_eco2 );
eco[3,n]~ordered_logistic(D_eco3 * lat_f[n,1] ,D_eco3 * cutoffs_eco3 );
rac[1,n]~bernoulli_logit( D_rac1 * lat_f[n,3] +D_rac1 * cutoffs_rac1 );
rac[2,n]~bernoulli_logit( D_rac2 * lat_f[n,3] +D_rac2 * cutoffs_rac2 );
rac[3,n]~bernoulli_logit( D_rac3 * lat_f[n,3] +D_rac3 * cutoffs_rac3 );
rac[4,n]~ordered_logistic(D_rac4 * lat_f[n,3] ,D_rac4 * cutoffs_rac4 );
rac[5,n]~ordered_logistic(D_rac5 * lat_f[n,3] ,D_rac5 * cutoffs_rac5 );
alt[1,n]~ordered_logistic(D_alt1 * lat_f[n,4] ,D_alt1 * cutoffs_alt1 );
alt[2,n]~ordered_logistic(D_alt2 * lat_f[n,4] ,D_alt2 * cutoffs_alt2 );
alt[3,n]~ordered_logistic(D_alt3 * lat_f[n,4] ,D_alt3 * cutoffs_alt3 );
imm[1,n]~ordered_logistic(D_imm1 * imm_lat[n] ,D_imm1 * cutoffs_imm1 );
imm[2,n]~ordered_logistic(D_imm2 * imm_lat[n] ,D_imm2 * cutoffs_imm2 );
imm[3,n]~ordered_logistic(D_imm3 * imm_lat[n] ,D_imm3 * cutoffs_imm3 );
imm[4,n]~ordered_logistic(D_imm4 * imm_lat[n] ,D_imm4 * cutoffs_imm4 );}
}

generated quantities{
cov_matrix[4] Omega;
real sd_eco;
real sd_com;
real sd_rac;
real sd_alt;
real rho_eco_com;
real rho_eco_rac;
real rho_eco_alt;
real rho_com_rac;
real rho_com_alt;
real rho_rac_alt;

Omega = quad_form_diag(L * L' ,sd_lat_f);

sd_eco=sqrt(Omega[1,1]);
sd_com=sqrt(Omega[2,2]);
sd_rac=sqrt(Omega[3,3]);
sd_alt=sqrt(Omega[4,4]);
rho_eco_com=Omega[2,1]/(sd_eco*sd_com);
rho_eco_rac=Omega[3,1]/(sd_eco*sd_rac);
rho_eco_alt=Omega[4,1]/(sd_eco*sd_alt);
rho_com_rac=Omega[3,2]/(sd_com*sd_rac);
rho_com_alt=Omega[4,2]/(sd_com*sd_alt);
rho_rac_alt=Omega[4,3]/(sd_rac*sd_alt);
}

