data {
  int<lower=0> N;
  int y [N];
  int x [N];
}

parameters {
  vector<lower=0> [max(x)] mu;
  vector<lower=0> [max(x)] sigma;
}

model {
  mu ~ normal(15, 5);
  sigma ~ normal(3, 2);
  y ~ neg_binomial_2(mu[x], exp(sigma[x]));
}
