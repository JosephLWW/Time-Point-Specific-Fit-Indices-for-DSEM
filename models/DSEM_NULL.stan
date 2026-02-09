data {
  int<lower=0> N;
  int<lower=0> Nt;
  matrix[N,6] y[Nt];
}

parameters {
  matrix[6,Nt]       ty;
  vector<lower=0>[6] sigmaeps; //resvar for obs var
}

model {

  for(i in 1:N){
    for(j in 1:Nt){
      for(k in 1:6){y[j,i,k] ~ normal(ty[k,j], sigmaeps[k]);}
    }
  }

  //priors
  //lu ~ normal(0,1);
  for(k in 1:6){ty[k,] ~ normal(0,10);}
  sigmaeps ~ normal(0, 1);

}


