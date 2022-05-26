// Complete pooling model

data {  				// Data block
  int<lower=0> N;		// Sample size, # of annual rates of change
  real dbio[N];	// Dependent variable, Annual rate of change
}

parameters {			// Parameter block
  real beta;		// Coefficient vector, pooled rate of change
  real<lower=0> sigma;	// Error scale
}
model {				// Model block
  dbio ~ normal(beta, sigma);
}
