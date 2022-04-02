# Bayesian Ordinal Quantile Regression with Bayesian Adaptive Lasso

BOQR-BAL (Bayesian Ordinal Quantile Regression with Bayesian Adaptive Lasso) performs quantile regression for ordinal outcomes with variable selection. Input data should be prepared as a matrix or dataframe with each row corresponding to an observation. The first column should be the ordinal outcome variable (coded with lowest value equal to 0), and the second column should be all 1s. The sampler can then be run for ```data``` at any quantile ```q``` between 0 and 1 with the following command:

```
result <- BOQR_BAL(data,q)
```

The output is the MCMC chain for the coefficients, with the first column corresponding to the intercept and the rest in order of the covariates as inputted via the original data. 
