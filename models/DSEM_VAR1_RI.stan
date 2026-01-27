data {
  int<lower=0> N;
  int<lower=0> Nt;
  matrix[N,6] y[Nt];
  //vector[Nt*6] x_mean;
  //matrix[Nt*6,Nt*6] x_cov;

}

parameters {
  vector<lower=0>[4] ly;
  vector[4]          ty;
  //vector<lower=0>[2] lu;
  vector<lower=0>[6] sigmaeps; //resvar for obs var
  vector<lower=0>[2] sigmazeta1; //resvar for obs var
  vector<lower=0>[2] sigmazeta2; //resvar for obs var
  vector[2]          ka; //intercepts
  matrix[2,2]        beta1; //ar(1) coefficients
  corr_matrix[2] phi1;
  corr_matrix[2] phi2;
  //real  eta[N,Nt,2];
  matrix[N,2] eta[Nt];
  
  // CAMBIO: Parametrizacion no centrada (z-scores)
  // Antes: matrix[N,2] u0;
  matrix[N,2] u0_z; 
}

transformed parameters{
  matrix[2,2] sigmaeta1;
  matrix[2,2] sigmaeta2;
  
  // CAMBIO: Recuperamos u0 real aqui
  matrix[N,2] u0;
  matrix[2,2] L_sigmaeta2;

  sigmaeta1 = quad_form_diag(phi1,sigmazeta1);
  sigmaeta2 = quad_form_diag(phi2,sigmazeta2);


  // Construccion no centrada
  L_sigmaeta2 = cholesky_decompose(sigmaeta2);
  // u0 = z * L' (transforma normales estandar a multivariadas correlacionadas)
  u0 = u0_z * L_sigmaeta2';
}

model {
  matrix[N,2] mueta[Nt];
  matrix[N,6] muy[Nt];
  // vector de mu0 no se usa explicitamente
  vector[2]   lu= rep_vector(1, 2);  
  
  // mean structure for eta, simple AR(1), no crossloadings
  for(i in 1:N){
    //t==1
    mueta[1,i,1] = ka[1]+lu[1]*u0[i,1];
    mueta[1,i,2] = ka[2]+lu[2]*u0[i,2];
  
    //t==rest
    for(j in 2:Nt){
      mueta[j,i,1] = ka[1]+lu[1]*u0[i,1]+beta1[1,1]*(eta[j-1,i,1]-lu[1]*u0[i,1]-ka[1])+beta1[1,2]*(eta[j-1,i,2]-lu[2]*u0[i,2]-ka[2]);
      mueta[j,i,2] = ka[2]+lu[2]*u0[i,2]+beta1[2,1]*(eta[j-1,i,1]-lu[1]*u0[i,1]-ka[1])+beta1[2,2]*(eta[j-1,i,2]-lu[2]*u0[i,2]-ka[2]);
    }
    
    for(j in 1:Nt){
      muy[j,i,1] = eta[j,i,1];
      muy[j,i,2] = ty[1]+ly[1]*eta[j,i,1];
      muy[j,i,3] = ty[2]+ly[2]*eta[j,i,1];
      
      muy[j,i,4] = eta[j,i,2];
      muy[j,i,5] = ty[3]+ly[3]*eta[j,i,2];
      muy[j,i,6] = ty[4]+ly[4]*eta[j,i,2];
      
    }
    
  }

  // Prior simple on scores z
  to_vector(u0_z) ~ std_normal();

  for (i in 1:N) {
    // sampling of u0 up
    
    for(j in 1:Nt){
      for(k in 1:6){y[j,i,k] ~ normal(muy[j,i,k], sigmaeps[k]);}

      eta[j,i,1:2] ~ multi_normal(mueta[j,i,1:2],sigmaeta1);    
    }
    
  }

  //priors
  //lu ~ normal(0,1);
  ty ~ normal(0,10);
  ly ~ normal(0,1);
  ka ~ normal(0,1);
  for(k in 1:2){beta1[,k] ~ normal(0,1);}
  sigmaeps ~ cauchy(0,5);
  sigmazeta1 ~ cauchy(0,5);
  sigmazeta2 ~ cauchy(0,5);
  phi1 ~ lkj_corr(1);
  phi2 ~ lkj_corr(1);
  

}


