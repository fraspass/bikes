# Title:        Utility functions for Hawkes processes (Jessica Ji's MSc project)
# Author:       Francesco Sanna Passino
# Affiliation:  Department of Mathematics, Imperial College London
# Email:        f.sannapassino@imperial.ac.uk
# Last updated: 17th July, 2021 - 10.00am

## Calculate intensity function
intensity = function(t,h,t_prime,params,check_t_prime=TRUE){
  ## Parameters
  lambda=params[1];  beta=params[2]; theta=params[3]; alpha=params[4]; delta=params[5]
  ## 
  if(check_t_prime==TRUE){
    t_prime = t_prime[t_prime < t]
  }
  return(lambda + sum(beta * exp(-theta * (t-h))) + sum(alpha * exp(-delta * (t-t_prime))))
}

## Simulate Hawkes process with scaled exponential self- and mutual-excitation function
hawkes_simulator = function(params, t_prime, T, lambda_star){
  ## Parameters
  lambda=params[1];  beta=params[2]; theta=params[3]; alpha=params[4]; delta=params[5]
  ## Initial points
  t_star = 0
  h = numeric(0)
  t_auxiliary = numeric(0)
  index = 0
  ## Simulate Hawkes process
  while(t_star < T){
    ## Propose new arrival time
    t_star = t_star - log(runif(1)) / lambda_star ## Important to pick a good upper bound lambda*
    if(t_star > T){
      break
    }
    if(index < length(t_prime)){
      if(t_star > t_prime[index+1]){
        for(j in (index+1):length(t_prime)){
          condition = (t_star < t_prime[j])
          if(condition == TRUE){
            index = j
            break
          } else {
            t_auxiliary = c(t_auxiliary, t_prime[j])
          }
        }
      }
    }
    lambda_prob = intensity(t_star,h,t_prime=t_auxiliary,
                            params=params,check_t_prime=FALSE) / lambda_star
    accept = sample(c(1,0), size=1, prob=c(lambda_prob, 1-lambda_prob))
    if(accept == 1){
      h = c(h, t_star)
    }
  }
  return(h)
}

## Transform parameters
transform_parameters = function(params, to_constrained=TRUE){
  tp = numeric(5)
  if(to_constrained==FALSE){
    tp[1] = log(params[1])
    tp[2] = log(params[2])
    tp[3] = log(params[3]-params[2])
    tp[4] = log(params[4])
    tp[5] = log(params[5]-params[4])
  } else {
    tp[1] = exp(params[1])
    tp[2] = exp(params[2])
    tp[3] = exp(params[2]) + exp(params[3])
    tp[4] = exp(params[4])
    tp[5] = exp(params[4]) + exp(params[5])
  }
  return(tp)
}

test = FALSE
if(test==TRUE){
  ## Parameters
  params = c(0.1, 0.25, 1.0, 0.1, 0.25)
  T = 10000
  ## Simulate t_prime as a Poisson process
  rho = 0.25
  set.seed(111)
  t_prime = rexp(n=1,rate=rho)
  while(TRUE){
    proposal = tail(t_prime,1)+rexp(n=1,rate=rho)
    if(proposal < T){
      t_prime = c(t_prime, proposal)
    } else {
      break 
    }
  }
  ## Simulate Hawkes process
  t = hawkes_simulator(params, t_prime, T, lambda_star=10)
}
