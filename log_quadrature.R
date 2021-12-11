library(tidyverse)
library(fitCensoRd)
library(microbenchmark)
library(data.table)

set.seed(242424241)



## Generates data
## Y|X ~ N(0,1)
## X ~ Weibull with parameters dependent on Bernoulli r.v.
## Censoring mechanism ~ Weibull dependent on Bernoulli r.v.
## Observed variable T is minimum of X and C

gen_dat <- function() {
N <- 1000
e <- rnorm(n = N, mean = 0, sd = 1)
z <- rbernoulli(n = N, p = 0.5) %>% as.numeric()
c1 <- rweibull(n = N, shape = 1, scale = 1.5)
c2 <- rweibull(n = N, shape = 2, scale = 1)
c <- ifelse(z==0, c1, c2)
x1 <- rweibull(n = N, shape = .75, scale = .25)
x2 <- rweibull(n = N, shape = 2, scale = .25)
x <- ifelse(z==0, x1, x2)
d <- as.numeric(x <= c)
sum(x <= c) ## Quantify proportion of censoring
t <- pmin(x, c)
y <- 1 + 0.5 * x + e

df <- matrix(c(y, x, c, t, d), ncol = 5)
colnames(df) <- c("y", "x", "c", "t", "d")
return(df)
}


## Y|X * X Likelihood, where X ~ Normal
## Method can be easily adapted to different densities for X
log_like <- function(df1, mean_y, df2, mean_x, sig_x, sig_y) {
  dnorm(x = df1, mean = mean_y, sd = sig_y) * dnorm(x = df2, mean = mean_x, sd = sig_x)
}

## MLE estimates under right-censored covariates
## Uses trapezoidal rule to estimate integral needed to calculate likelihood

loglik_quad <- function(params, df) {
  
  # data.table is used for speed of group sums (compared to something like tapply)
  require(data.table)
  
  # Parameters to be estimated
  beta0 <- params[1]
  beta1 <- params[2]
  sig_y <- params[3]
  
  
  # Censored covariate parameters
  mu_x <- params[4]
  sig_x <- params[5]
  
  ## Observed-data likelihood
  
  obs <- df[, "x"] <= df[, "c"] # Create matrix of observed data
  mu_y <- beta0 + beta1 * df[obs, "x"]
  obs_like <- log_like(df[obs, "y"], mu_y, df[obs, "x"], mu_x, sig_x, sig_y)
  ll_obs <- sum(log(obs_like))
  
  ## Censored likelihood
  
  range_t <- max(df[obs,"x"]) - df[!obs, "t"] # Calculates support of each integral for each unobserved i
  
  inelig <- range_t < 0 # Checks that greatest observed X is larger than censored T_i
  range_t <- range_t[!inelig]
  n_cens <- length(df[!obs, "t"][!inelig]) # Number of (eligible) censored observations
  
  split <- rep.int(df[!obs, "t"][!inelig], rep.int(101, n_cens)) + rep.int(seq.int(from = 0, to = 1, by = .01), n_cens) * 
    rep.int(range_t, rep.int(101, n_cens)) # Partitions each individual support into 101 equal partitions
                                           # 101 equal partitions might seem odd but want to avoid floating point errors and dividing by 100 seems easier
  
  expand_df <- cbind(split, rep.int(df[!obs, "y"][!inelig], rep.int(101, n_cens)), # Expands out all pairwise comparisons of each partitioned x
                     rep.int(1:n_cens, rep.int(101, n_cens)))                      # for each Y_i
  
  vals <- log_like(expand_df[, 2], beta0 + beta1*expand_df[, 1], expand_df[, 1], mu_x, sig_x, sig_y) * 
    rep.int(c(0.5, rep(1, 99), 0.5), times = n_cens) * rep.int(range_t/100, rep.int(101, n_cens))
    # Multiply by 0.5 for the 1st and last entry as per trapezoidal rule, as we used equal partitions
  
  new_df <- cbind(vals, expand_df[, 3]) # Combine results of previous and grouping index in preparation for grouped sum
  dt <- data.table(new_df) # data.table used here for speed
  setkey(dt,V2) # used by data.table
  dt <- dt[, sum(vals), by = V2] 
  ll_censor <- sum(log(dt$V1))
  ll <- ll_obs + ll_censor # Final log-likelihood
  return(-ll)
}

## Uses nlm to minimize log-likelihood calculated in previous
linreg_quad <- function(params0, data, steptol = 1e-6, iterlim = 100) {

  suppressWarnings(

    ols <- nlm(f = loglik_quad, p = params0, steptol = steptol, iterlim = iterlim, hessian = TRUE, df = data)
)
  
  param_est <- ols$estimate
  param_se <- c(sqrt(diag(solve(ols$hessian)))[1:2], NA, sqrt(diag(solve(ols$hessian)))[4], NA)
  param_df <- data.frame(est = param_est, se = param_se)
  rownames(param_df) <- c("beta0", "beta1", "sigmaY", "alpha0", "sigmaX")
  return(list(coeff = param_df, code = ols$code))
}



## Compares estimates to Atem & Matsouaka
estimates <- list()

for(i in 1:1000) {
  set.seed(i)
  df <- gen_dat()
  params = c(0, 0, 1, 0, 1)
  fit1 <- linreg_quad(params0 = params, data = df) %>% data.frame()
  fit2 <- fitCensoRd::linreg_am(params0 = params, Y = "y", X = "x", C = "c", data = df) %>% data.frame()
  fit1$seed <- i
  fit2$seed <- i
  results <- data.frame(fit1, fit2)
  estimates[[i]] <- results
}

estimates_list = do.call(rbind, estimates)

## Measure (absolute) difference in parameter estimates/standard errors
estimates_list$diffest <- abs(estimates_list$coeff.est - estimates_list$coeff.est.1)
estimates_list$diffse <- estimates_list$coeff.se - estimates_list$coeff.se.1


write_csv(estimates_list, file = "")



