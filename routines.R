# Quantile functions
ql <- function(p) quantile(p, 0.025)
qu <- function(p) quantile(p, 0.975)


######################################################################################
######################################################################################
######################################################################################
# Weibull
######################################################################################
######################################################################################
######################################################################################


#--------------------------------------------------------------------------------------------------
# Weibull -log-posterior function
#--------------------------------------------------------------------------------------------------

log_postW <- function(par){
  sigma <- exp(par[1]); nu <- exp(par[2]);
  
  # Terms in the log log likelihood function
  ll_haz <- sum(hweibull(t_obs, sigma, nu, log = TRUE))
  
  ll_chaz <- sum(chweibull(survtimes, sigma, nu))
  
  log_lik <- -ll_haz + ll_chaz 
  
  # Log prior
  
  log_prior <- -dgamma(sigma, shape = 2, scale = 2, log = TRUE) - 
    dgamma(nu, shape = 2, scale = 2, log = TRUE) 
  
  # Log Jacobian
  
  log_jacobian <- -par[1] - par[2] 
  
  # log posterior
  
  log_post0 <- log_lik + log_prior + log_jacobian
  
  return(as.numeric(log_post0))
}

######################################################################################
######################################################################################
######################################################################################
# Power Generalised Weibull
######################################################################################
######################################################################################
######################################################################################


#--------------------------------------------------------------------------------------------------
# Power Generalised Weibull -log-posterior function
#--------------------------------------------------------------------------------------------------

log_postPGW <- function(par){
  sigma <- exp(par[1]); nu <- exp(par[2]); gamma <- exp(par[3])
  
  # Terms in the log log likelihood function
  ll_haz <- sum(hpgw(t_obs, sigma, nu, gamma, log = TRUE))
  
  ll_chaz <- sum(chpgw(survtimes, sigma, nu, gamma))
  
  log_lik <- -ll_haz + ll_chaz 
  
  # Log prior
  
  log_prior <- -dgamma(sigma, shape = 2, scale = 2, log = TRUE) - 
    dgamma(nu, shape = 2, scale = 2, log = TRUE) - 
    dgamma(gamma, shape = 2, scale = 2, log = TRUE) 
  
  # Log Jacobian
  
  log_jacobian <- -par[1] - par[2] - par[3]
  
  # log posterior
  
  log_post0 <- log_lik + log_prior + log_jacobian
  
  return(as.numeric(log_post0))
}

######################################################################################
######################################################################################
######################################################################################
# Logistic ODE
######################################################################################
######################################################################################
######################################################################################

#-----------------------------------------------------------------------------
# Logistic ODE Hazard Function: Analytic solution
#-----------------------------------------------------------------------------
# t: time (positive)
# lambda: intrinsic growth rate (positive)
# kappa: upper bound (positive)
# h0: hazard initial value (positive)
hlogisode <- function(t, lambda, kappa, h0, log = FALSE){
  lhaz <-  log(kappa) + log(h0) - log( (kappa-h0)*exp(-lambda*t) + h0)
  if (log)  return(lhaz)
  else return(exp(lhaz))
}

#-----------------------------------------------------------------------------
# Logistic ODE Cumulative Hazard Function: Analytic solution
#-----------------------------------------------------------------------------
# t: time (positive)
# lambda: intrinsic growth rate (positive)
# kappa: upper bound (positive)
# h0: hazard initial value (positive)
chlogisode <- function(t, lambda, kappa, h0){
  chaz <-   kappa*( log( (kappa-h0)*exp(-lambda*t) + h0 ) - log(kappa) + lambda*t )/lambda
  return(chaz)
}

#-----------------------------------------------------------------------------
# Logistic ODE random number generation
#-----------------------------------------------------------------------------
# t: time (positive)
# lambda: intrinsic growth rate (positive)
# kappa: upper bound (positive)
# h0: hazard initial value (positive)
rlogisode <- function(n, lambda, kappa, h0){
  u <- runif(n)
  times <- log( 1 + kappa*( exp(-lambda*log(1-u)/kappa) -1 )/h0    )/lambda 
  return(as.vector(times))
}

#-----------------------------------------------------------------------------
# Logistic ODE Probability Density Function: Analytic solution
#-----------------------------------------------------------------------------
# t: time (positive)
# lambda: intrinsic growth rate (positive)
# kappa: upper bound (positive)
# h0: hazard initial value (positive)
dlogisode <- function(t, lambda, kappa, h0, log = FALSE){
  lden <-  hlogisode(t, lambda, kappa, h0, log = TRUE) - 
    chlogisode(t, lambda, kappa, h0)
  if (log)  return(lden)
  else return(exp(lden))
}



#--------------------------------------------------------------------------------------------------
# Logistic ODE -log-likelihood function: Analytic solution
#--------------------------------------------------------------------------------------------------

log_likL <- function(par){
  lambda <- exp(par[1]); kappa <- exp(par[2]); h0 <- exp(par[3])
    # Terms in the log log likelihood function
    ll_haz <- sum(hlogisode(t_obs, lambda, kappa, h0, log = TRUE))
    
    ll_chaz <- sum(chlogisode(survtimes, lambda, kappa, h0 ))
    
    log_lik <- -ll_haz + ll_chaz 
    
    return(log_lik)
}


#--------------------------------------------------------------------------------------------------
# Logistic ODE -log-posterior function: Analytic solution
#--------------------------------------------------------------------------------------------------

log_postL <- function(par){
  lambda <- exp(par[1]); kappa <- exp(par[2]); h0 <- exp(par[3])
  
  # Terms in the log log likelihood function
  ll_haz <- sum(hlogisode(t_obs, lambda, kappa, h0, log = TRUE))
  
  ll_chaz <- sum(chlogisode(survtimes, lambda, kappa, h0 ))
  
  log_lik <- -ll_haz + ll_chaz 
  
  # Log prior
  
  log_prior <- -dgamma(exp(par[1]), shape = 2, scale = 2, log = TRUE) - 
    dgamma(exp(par[2]), shape = 2, scale = 2, log = TRUE) -
    dgamma(exp(par[3]), shape = 2, scale = 2, log = TRUE) 
  
  # Log-Jacobian
  
  log_jacobian <- - par[1] - par[3] - par[3]
    
  # log posterior
  
  log_post0 <- log_lik + log_prior + log_jacobian

  return(as.numeric(log_post0))
}



######################################################################################
######################################################################################
######################################################################################
# Hazard-Response 
######################################################################################
######################################################################################
######################################################################################



#--------------------------------------------------------------------------------------------------
# Hazard-Response ODE Solver 
#--------------------------------------------------------------------------------------------------


# Hazard-Response ODE
hazmodHR <- function(t, y, par) {
  # state variables
  h <- y[1]
  q <- y[2]
  CH <- y[3]

  # parameters
  lambda <- par[1]
  kappa <- par[2]
  alpha <- par[3]
  beta <- par[4]

  # model equations
  dh <-  lambda*h*(1 - h/kappa) - alpha*q*h
  dq <-  beta*q*(1-q/kappa) - alpha*q*h
  dCH <- h

  # result
  return( list(c(dh, dq, dCH)) )

}

#--------------------------------------------------------------------------------------------------
# Hazard-Response ODE Solver : log h and log q
#--------------------------------------------------------------------------------------------------


# Hazard-Response ODE
hazmodHRL <- function(t, y, par) {
  # state variables
  lh <- y[1]
  lq <- y[2]
  CH <- y[3]
  
  # parameters
  lambda <- par[1]
  kappa <- par[2]
  alpha <- par[3]
  beta <- par[4]
  
  # model equations
  dlh <-  lambda*(1 - exp(lh)/kappa) - alpha*exp(lq)
  dlq <-  beta*(1-exp(lq)/kappa) - alpha*exp(lh)
  dCH <- exp(lh)
  
  # result
  return( list(c(dlh, dlq, dCH)) )
  
}

#--------------------------------------------------------------------------------------------------
# Hazard-Response ODE : Jacobian of the system of ODEs
#--------------------------------------------------------------------------------------------------

jacODE <- function(t, y, par) {
  # state variables
  h <- y[1]
  q <- y[2]
  CH <- y[3]
  
  # parameters
  lambda <- par[1]
  kappa <- par[2]
  alpha <- par[3]
  beta <- par[4]
  
  J1 <- c(lambda*(1-2*h/kappa) - alpha*q, -alpha*h, 0)
  J2 <- c(-alpha*q, beta*(1-2*q/kappa) - alpha*h, 0)
  J3 <- c(1, 0, 0)
  
  jacob <- as.matrix(rbind(J1, J2, J3))
  return(jacob)
}



#--------------------------------------------------------------------------------------------------
# Hazard-Response ODE : Jacobian of the system of ODEs
#--------------------------------------------------------------------------------------------------

jacODEL <- function(t, y, par) {
  # state variables
  lh <- y[1]
  lq <- y[2]
  CH <- y[3]
  
  # parameters
  lambda <- par[1]
  kappa <- par[2]
  alpha <- par[3]
  beta <- par[4]
  
  J1 <- c( - lambda*exp(lh)/kappa, - alpha*exp(lq), 0)
  J2 <- c(- alpha*exp(lh),- beta*exp(lq)/kappa, 0)
  J3 <- c(exp(lh), 0, 0)
  
  jacob <- as.matrix(rbind(J1, J2, J3))
  return(jacob)
}

#--------------------------------------------------------------------------------------------------
# Hazard-Response ODE -log-likelihood function: Solver 
#--------------------------------------------------------------------------------------------------

log_likHR <- function(par){
  # Numerical solution for the ODE
  params  <- c(lambda = exp(par[1]), kappa = exp(par[2]), alpha = exp(par[3]),
               beta = exp(par[4]))
  init      <- c(h = h0, q = q0, H = 0 )
  times  <- c(0,survtimes) 
  out <- ode(init, times, hazmodHR, params, method = "lsode", 
             jacfunc = jacODE, jactype = "fullusr")[-1,]
  
  # Terms in the log log likelihood function
  ll_haz <- sum(log(as.vector(out[status,2])))
  
  ll_chaz <- sum(as.vector(out[,4]))
  
  log_lik <- -ll_haz + ll_chaz
  
  return(log_lik)
}



#--------------------------------------------------------------------------------------------------
# Hazard-Response ODE -log-posterior function: Solver
#--------------------------------------------------------------------------------------------------

log_postHR <- function(par){
  # Numerical solution for the ODE
  params  <- c(lambda = exp(par[1]), kappa = exp(par[2]), alpha = exp(par[3]),
               beta = exp(par[4]))
  init      <- c(h = h0, q = q0, H = 0 )
  times  <- c(0,survtimes) 
  out <- ode(init, times, hazmodHR, params, method = "lsode", 
             jacfunc = jacODE, jactype = "fullusr")[-1,]
  
  if(any(out[,2]<=0)){ return(1e6) }
  
  else{
    # Terms in the log log likelihood function
    ll_haz <- sum(log(as.vector(out[status,2])))
    
    ll_chaz <- sum(as.vector(out[,4]))
    
    log_lik <- -ll_haz + ll_chaz
    
    # Log prior
    
    log_prior <- -dgamma(exp(par[1]), shape = 2, scale = 2, log = TRUE) -
      dgamma(exp(par[2]), shape = 2, scale = 2, log = TRUE) -
      dgamma(exp(par[3]), shape = 2, scale = 2, log = TRUE) -
      dgamma(exp(par[4]), shape = 2, scale = 2, log = TRUE) 
    
    # Log Jacobian
    
    log_jacobian <- - par[1] - par[2] - par[3] - par[4] 
    
    # log posterior
    
    log_post0 <- log_lik + log_prior + log_jacobian
    
    return(as.numeric(log_post0))
  }
}


#--------------------------------------------------------------------------------------------------
# Hazard-Response ODE -log-likelihood function: Solver (log h and log q)
#--------------------------------------------------------------------------------------------------

log_likHRL <- function(par){
  # Numerical solution for the ODE
  params  <- c(lambda = exp(par[1]), kappa = exp(par[2]), alpha = exp(par[3]),
               beta = exp(par[4]))
  init      <- c(lh = log(h0), lq = log(q0), H = 0 )
  times  <- c(0,survtimes) 
  out <- ode(init, times, hazmodHRL, params, method = "lsode", 
             jacfunc = jacODEL, jactype = "fullusr")[-1,]
  
  # Terms in the log log likelihood function
  ll_haz <- sum(as.vector(out[status,2]))
  
  ll_chaz <- sum(as.vector(out[,4]))
  
  log_lik <- -ll_haz + ll_chaz
  
  return(log_lik)
}



#--------------------------------------------------------------------------------------------------
# Hazard-Response ODE -log-posterior function: Solver (log h and log q)
#--------------------------------------------------------------------------------------------------

log_postHRL <- function(par){
  # Numerical solution for the ODE
  params  <- c(lambda = exp(par[1]), kappa = exp(par[2]), alpha = exp(par[3]),
               beta = exp(par[4]))
  init      <- c(lh = log(h0), lq = log(q0), H = 0 )
  times  <- c(0,survtimes) 
  out <- ode(init, times, hazmodHRL, params, method = "lsode", 
             jacfunc = jacODEL, jactype = "fullusr")[-1,]
  

    # Terms in the log log likelihood function
    ll_haz <- sum(as.vector(out[status,2]))
    
    ll_chaz <- sum(as.vector(out[,4]))
    
    log_lik <- -ll_haz + ll_chaz
    
    # Log prior
    
    log_prior <- -dgamma(exp(par[1]), shape = 2, scale = 2, log = TRUE) -
      dgamma(exp(par[2]), shape = 2, scale = 2, log = TRUE) -
      dgamma(exp(par[3]), shape = 2, scale = 2, log = TRUE) -
      dgamma(exp(par[4]), shape = 2, scale = 2, log = TRUE) 
    
    # Log Jacobian
    
    log_jacobian <- - par[1] - par[2] - par[3] - par[4] 
    
    # log posterior
    
    log_post0 <- log_lik + log_prior + log_jacobian
    
    return(as.numeric(log_post0))
}



#--------------------------------------------------------------------------------------------------
# Hazard-Response ODE -log-posterior function: Solver. Unknown initial conditions
#--------------------------------------------------------------------------------------------------

log_postHRX <- function(par){
  # Numerical solution for the ODE
  params  <- c(lambda = exp(par[1]), kappa = exp(par[2]), alpha = exp(par[3]),
               beta = exp(par[4]), h0 = exp(par[5]), q0 = exp(par[6]))
  init      <- c(h = h0, q = q0, H = 0 )
  times  <- c(0,survtimes) 
  out <- ode(init, times, hazmodHR, params, method = "lsode", 
             jacfunc = jacODE, jactype = "fullusr")[-1,]
  
  if(any(out[,2]<=0)){ return(1e6) }
  
  else{
    # Terms in the log log likelihood function
    ll_haz <- sum(log(as.vector(out[status,2])))
    
    ll_chaz <- sum(as.vector(out[,4]))
    
    log_lik <- -ll_haz + ll_chaz
    
    # Log prior
    
    log_prior <- -dgamma(exp(par[1]), shape = 2, scale = 2, log = TRUE) -
      dgamma(exp(par[2]), shape = 2, scale = 2, log = TRUE) -
      dgamma(exp(par[3]), shape = 2, scale = 2, log = TRUE) -
      dgamma(exp(par[4]), shape = 2, scale = 2, log = TRUE) -
      dgamma(exp(par[5]), shape = 2, scale = 5e-3, log = TRUE) -
      dgamma(exp(par[6]), shape = 2, scale = 5e-7, log = TRUE) 
      
    
    # Log Jacobian
    
    log_jacobian <- - par[1] - par[2] - par[3] - par[4] - par[5] - par[6]
    
    # log posterior
    
    log_post0 <- log_lik + log_prior + log_jacobian
    
    return(as.numeric(log_post0))
  }
}

#--------------------------------------------------------------------------------------------------
# Random number generation (Hazard-Response ODE)
# Use gamma = kappa and delta = alpha to obtain the model in the manuscript
#--------------------------------------------------------------------------------------------------
# n : number of simulations
# par: parameters (lambda, kappa, alpha, beta, h0, q0)
# tmin: minimum time to search for a solution for the simulation (typically tmin = 0)
# tmax: maximum time to search for a solution for the simulation
# step: step between nodes (equidistant)
# Returns the simulated times, the interpolating function, and the ODESolver object

simHR <- function(n, par, tmin, tmax, step){
  # Parameters
  params  <- c(lambda = par[1], kappa = par[2], alpha = par[3],
               beta = par[4])
  init0      <- c(h = h0, q = q0, H = 0 )
  # y-axis coordinates
  times = seq(tmin, tmax, by = step)
  # ODE solution
  sol <- ode(init0, times, hazmodHR, params, method = "lsode")
  # Interpolating function for the inverse of the cumulative hazard
  invH <- cm.splinefun(x = sol[,4], y = sol[,1])
  invHv <- Vectorize(function(t) invH(t))
  
  sim <- invHv(-log(runif(n)))
  out <- list( sim = sim, ODEsol = sol, invH = invHv )
  return(out)
}
