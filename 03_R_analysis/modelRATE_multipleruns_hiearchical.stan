data {
  int<lower=1> N; 
  int<lower=1> N_group;   // Number of modelruns
  int<lower=1,upper=N_group> modelruns[N];
  vector[N] dbio;
} 
parameters {
  vector[N_group] a1;
  real mu_a1;
  real<lower=0> sigma_a1;
  real<lower=0> sigma_y;
}
model {
  vector[N] y_hat;
  for (i in 1:N)
    y_hat[i] = a1[modelruns[i]];

  mu_a1 ~ normal(0, 1);
  sigma_a1 ~ cauchy(0, 1);
  sigma_y ~ cauchy(0, 1);
  a1 ~ normal(mu_a1, sigma_a1);
  dbio ~ normal(y_hat, sigma_y);
}
