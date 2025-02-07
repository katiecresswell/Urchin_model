//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//


functions {
  real [ ] findbiomass(real a, real b, real c, real [ ] weight, real [ , ] Gi, real [ ] Histcatch, real futurecatch, int [ ] xt, int finalt, real reefa, int abc) {
    real  Totalbiomass[finalt];
    real  Detectbiomass[finalt];
    real  Detectdensity[finalt];
    real  Exploitablebiomass[finalt];
    real  sumcatch[finalt];
    real  Nsize[50,finalt];
    real  Ncatch[50,finalt];
    real  Nrecruit[finalt];
    int   j[50];
    real  catchselectivity[50];
    real  recruitment;
    real  thiscatch;
    real  fcatch;
    int   catchmin;
    int   endt = max(xt);
    
    for (i in 1:50) {
      for (t in 1:endt) Ncatch[i,t] = 0;
      Nsize[i,1] = 0;
      // size min for catch
      catchmin = 24;
      if (i<catchmin) catchselectivity[i] = 0;
      else catchselectivity[i] = 1;
    }
    Detectbiomass[1] = 0.1;
    Totalbiomass[1] = 0.1;
    fcatch = .02;
    Detectdensity[1] = Detectbiomass[1]/reefa;
    Exploitablebiomass[1] = Detectdensity[1];
    sumcatch[1] = 0;
    Nrecruit[1] = c/(1+exp(-a*(1-b)));
    for (t in 2:endt){
      Detectbiomass[t] = 0;
      Exploitablebiomass[t] = 0;
      Totalbiomass[t] = 0;
      sumcatch[t] = 0;
      
      for (i in 1:50){
        // add recruits to first size only
        if (i==1) {
          recruitment = c/(1+exp(-a*(t-b)));
          Nrecruit[t] = recruitment;
          Nsize[i,t] = Gi[i,i]*(Nsize[i,t-1]);
        } else {
          recruitment = 0;
          // number size i that stay at size i
          Nsize[i,t] = Gi[i,i]*(Nsize[i,t-1]);
          for (jint in 1:(i-1)){
            //number at size <i that grow into size i
            Nsize[i,t] += Gi[i,jint]*Nsize[jint,t-1];
          }
        }
        // new num at size i * mortality then + recruits
        Nsize[i,t] = Nsize[i,t]*exp(-0.11) + recruitment;
        if (Nsize[i,t]<0) Nsize[i,t]=0;
        Exploitablebiomass[t] += (Nsize[i,t]*weight[i]*catchselectivity[i]);
      } 
      if (Exploitablebiomass[t]<0) Exploitablebiomass[t] = 0;
      thiscatch = 0;
      // if "never fished" scenario then all catches zero
      if (abc==7) {
        for (i in 1:50) Ncatch[i,t] = 0;
      }
      else {
        if (t>48){ // 49 is the year historical fishing began 2009
          if (t<=63) { // 63 is final year of historical catch 2023
            thiscatch = Histcatch[t-48];
            // if the catch (historical or future scenario) is less than exploitable biomass available then good 
            if (thiscatch<=Exploitablebiomass[t]) {
              if (Exploitablebiomass[t]>0) for (i in 1:50) Ncatch[i,t] = Nsize[i,t]*catchselectivity[i]*thiscatch/Exploitablebiomass[t];
              else for (i in 1:50) Ncatch[i,t] = 0;
            } else { // otherwise make it equal to exploitable biomass available 
              if (Exploitablebiomass[t]>0) for (i in 1:50) Ncatch[i,t] = Nsize[i,t]*catchselectivity[i];
              else for (i in 1:50) Ncatch[i,t] = 0;
            } 
          } else { // else if time > 2023 use future scenario for catch fcatch 0.02*expB
            thiscatch = (futurecatch-1)*fcatch*Exploitablebiomass[t];
            // error checking
            if (thiscatch>Exploitablebiomass[t]) print(" how ??? t ",t," exp ",Exploitablebiomass[t]," thiscatch ",thiscatch," futurec ",futurecatch);
            for (i in 1:50) Ncatch[i,t] = Nsize[i,t]*catchselectivity[i]*thiscatch/Exploitablebiomass[t];
          } // end t>2023
        }
      }
      // add up sizes into total and detectable biomass
      for (i in 1:50){    
        Nsize[i,t] = Nsize[i,t] - Ncatch[i,t];
        if (i>=10) Detectbiomass[t] += Nsize[i,t]*weight[i];
        Totalbiomass[t] += Nsize[i,t]*weight[i];
        sumcatch[t] += Ncatch[i,t]*weight[i];
      }
      // convert to density by diving by reef area for that region z
      Detectdensity[t] = Detectbiomass[t]/reefa;
    }
    for (t in 1:endt) {
      if (Detectdensity[t]<0) Detectdensity[t]=0;
      if (Exploitablebiomass[t]<0) Exploitablebiomass[t]=0;
      if (sumcatch[t]<0) sumcatch[t]=0;
    }
   
    if (abc==2) return (Exploitablebiomass[xt]);
    else if ((abc==3)||(abc==5)) return (sumcatch[xt]);
    else return (Detectdensity[xt]);
  } 
}
data {
  int<lower=0> N;  
  int<lower=0> fcmax; 
  int<lower=0> x[N]; 
  real<lower=0> Y[N]; 
  real<lower=0> G[50,50];
  real<lower=0> histcatch[15];
  real<lower=0> weight[50];
  real<lower=0> reefarea;
  //real<lower=0> cpram;
  
  // Variables for predictions
  int<lower=0> newt; // number of years of prediction data set
  int<lower=0> newx[newt]; // number of years of prediction data set
  
} 
parameters {
  real<lower=0> apram; //shape was 100
  real<lower=0> bpram;  // inflection
  real<lower=0> cpram;   //scale was 3,000,000
  real<lower=0> tau;
  
}  

transformed parameters {
  real<lower=0> m[N];
  real<lower=0> sigma;
  m = findbiomass(apram,bpram,cpram,weight,G, histcatch, 1, x, newt, reefarea, 1);
  sigma = 1/sqrt(tau);
}
model {
  real disp_sigma[N];
  // priors
  apram ~ normal(10,10); 
  bpram ~ normal(30,10); 
  cpram ~ normal(2000000,1000000); 
  tau ~ gamma(.0001,.0001);
  
  //likelihood
  for (n in 1:N) {
    Y[n] ~ normal(m[n],sigma) T[0,];
  }
}
generated quantities{
  real Y_meanfished[fcmax,newt]; 
  real Exploitbio[fcmax,newt];
  real Timecatch[fcmax,newt];
  real Y_predfished[fcmax,newt]; 
  real Y_neverfished[newt]; 
  
  // Posterior parameter distribution of the mean
  for (fc in 1:fcmax){
    Y_meanfished[fc,] = findbiomass(apram,bpram,cpram,weight,G, histcatch, fc, newx, newt, reefarea, 1);
    Exploitbio[fc,] = findbiomass(apram,bpram,cpram,weight,G, histcatch, fc, newx, newt, reefarea, 2);
    Timecatch[fc,] = findbiomass(apram,bpram,cpram,weight,G, histcatch, fc, newx, newt, reefarea, 3);
    if (fc==1) {
      Y_neverfished = findbiomass(apram,bpram,cpram,weight,G, histcatch, fc, newx, newt, reefarea, 7);
    }
    for(t in 1:newt){
      // Posterior predictive distribution
      Y_predfished[fc,t] = normal_rng(Y_meanfished[fc,t], sigma) ;
    }
  }
}


