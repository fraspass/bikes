# Title:        Utility functions for Hawkes processes
# Author:       Francesco Sanna Passino
# Affiliation:  Department of Mathematics, Imperial College London
# Email:        f.sannapassino@imperial.ac.uk
# Last updated: 17th July, 2021 - 8.50am

## Calculate intensity function
intensity = function(t,h,lambda,beta,theta){
  return(lambda + sum(beta * exp(-theta * (t-h))))
}

## Simulate Hawkes process
hawkes_simulator = function(lambda, beta, theta, T){
  ## Initial points
  t_star = 0
  h = numeric(0)
  ## Simulate Hawkes process
  while(t_star < T){
    lambda_star = intensity(t_star, h, lambda, beta, theta)
    ## Propose new arrival time
    t_star = t_star - log(runif(1)) / lambda_star
    if(t_star > T){
      break
    }
    lambda_prob = intensity(t_star,h,lambda, beta, theta) / lambda_star
    accept = sample(c(1,0), size=1, prob=c(lambda_prob, 1-lambda_prob))
    if(accept == 1){
      h = c(h, t_star)
    }
  }
  return(h)
}

## Transform parameters
transform_parameters = function(params,inverse=FALSE){
  tp = numeric(3)
  if(inverse==TRUE){
    tp[1] = log(params[1])
    tp[2] = log(params[2])
    tp[3] = log(params[3]-params[2])
  } else {
    tp[1] = exp(params[1])
    tp[2] = exp(params[2])
    tp[3] = exp(params[2]) + exp(params[3])
  }
  return(tp)
}

## Log-likelihood
log_likelihood = function(params,t,T){
  ## Reparametrise
  trans_params = transform_parameters(params, inverse=FALSE)
  lambda = trans_params[1]
  beta = trans_params[2]
  theta = trans_params[3]
  ## Calculate likelihood 
  logl = -lambda*T + beta/theta * sum(exp(-theta * (T-t)) - 1) 
  A = 0
  t = c(t,T)
  for(i in 1:length(t)){
    logl = logl + log(lambda + beta * A)
    A = exp(-theta*(t[i+1]-t[i])) * (1 + A)
  }
  return(logl)
}

## Optimiser
optim_logl = function(starting_values,t,T){
  neg_logl = function(starting_values) -log_likelihood(starting_values, t=t, T=T)
  return(optim(par=starting_values, fn=neg_logl))
}

## Compensator
compensator = function(t,h,lambda,beta,theta){
  return(lambda*t - beta / theta * sum(exp(-theta * (t-h[h<t])) - 1))
}

###############
### Testing ###
###############

run_test = FALSE
if(run_test == TRUE){
  ## Parameters
  lambda = 0.01
  beta = 0.1
  theta = 1
  T = 1000000
  ## Simulate
  set.seed(111)
  t = hawkes_simulator(lambda,beta,theta,T) ## Should generate 11117 events in ~10 seconds
  ## Maximum likelihood estimation
  res = optim_logl(starting_values=transform_parameters(c(lambda,beta,theta),inverse=TRUE),t=t,T=T)
  transform_parameters(res$par)
  ## Test for uniformity of p-values
  ## Warning: inefficient calculation (~2 seconds for 10k events)
  pv = 1-exp(-diff(sapply(t,function(x) compensator(t=x,h=t,lambda=lambda,beta=beta,theta=theta))))
  hist(pv,freq=FALSE,main='Histogram of p-values')
  qqplot(pv, seq(0,length(t)-1)/(length(t)-1),type='l',lwd=2,col='red',
         ylab='Theoretical quantiles',xlab='Observed quantiles', main='Q-Q plot')
  abline(0,1,lty='dotted')
  ks.test(pv,'punif') ## Test for uniformity: fail to reject null H0 --> good model fit
}
