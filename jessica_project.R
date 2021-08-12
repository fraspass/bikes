# Title:        Hawkes process for starting times of bike journeys
# Author:       Francesco Sanna Passino
# Affiliation:  Department of Mathematics, Imperial College London
# Email:        f.sannapassino@imperial.ac.uk

## Empty list object
weekly_data = list()
weekly_data_test = list()
## List all the .csv file names (and sort them in case they are not ordered by number)
file_names = sort(list.files('data/', pattern='2021.csv'))
## Total number of files
n_weeks = length(file_names)
## Import 4 weeks of data as training set
for(week in (n_weeks-5):(n_weeks-2)){
  weekly_data[[week]] = read.table(paste('data/',file_names[week],sep=''), 
                                   sep=',', header=FALSE, 
                                   col.names=c('start_id','end_id','start_time','duration'))
}
## Last 2 weeks as test set
for(week in (n_weeks-1):n_weeks){
  weekly_data_test[[week]] = read.table(paste('data/',file_names[week],sep=''), 
                                        sep=',', header=FALSE, 
                                        col.names=c('start_id','end_id','start_time','duration'))
}
## Create a unique dataframe
library(dplyr)
df = dplyr::bind_rows(weekly_data)
df = transform(df, end_time = start_time + duration)
## Import stations
stations = read.table('data/santander_locations.csv', sep=',', header=TRUE)
## Find StationID for Hyde Park Corner
id = stations[grepl('Hyde Park Corner', stations$StationName),]$Station.Id
## Filter 
t1 = df$start_time[df$start_id == id]
t1_prime = df$end_time[df$end_id == id]
## Obtain the start of the observation period
t0 = floor(min(c(min(t1),min(t1_prime))) / 60 / 60 / 24) * (60 * 60 * 24)
## Shift the times for convenience + add random uniform noise (and sort) 
## + transform times to hours (for convenience)
t = sort((t1 - t0) + runif(n=length(t1))) / 60 / 60
t_prime = sort(t1_prime - t0 + runif(n=length(t1_prime))) / 60 / 60
## Total time of observation -- 28 days * 24 hours = 672 hours (as expected -- 4 weeks of data)
T = ceiling(max(c(max(t),max(t_prime))) / 24) * 24

## *Negative* log-lihelihood
log_likelihood = function(params, t, t_prime, T, B_sequence){
  ## Reparametrise
  trans_params = transform_parameters(params, to_constrained=TRUE)
  lambda = trans_params[1]
  beta = trans_params[2]
  theta = trans_params[3]
  alpha = trans_params[4]
  delta = trans_params[5]
  ## B_function is a list with len(t) elements 
  ## B_sequence[[i]] contains the difference between t[i] and the t_primes such that:
  ##    t_prime < t[i] AND t_prime >= t[i-1]
  ## If B_sequence is not passed to the function, calculate it. 
  ## For optimisation, a B_sequence ** must ** be passed to the function for efficient optimisation. 
  if(missing(B_sequence) == TRUE){
    B_sequence = B_sequence_calculator(t, t_prime)
  }
  ## Calculate likelihood 
  logl = -lambda*T + beta/theta * sum(exp(-theta * (T-t)) - 1)
  logl = logl + alpha/delta * sum(exp(-delta * (T-t_prime)) - 1)
  A = 0
  B = sum(exp(-delta*B_sequence[[1]])) ## IMPORTANT: initial value of B is *not* 0
  for(i in 1:length(t)){
    logl = logl + log(lambda + beta * A + alpha * B)
    if(i < length(t)){
      A = exp(-theta*(t[i+1]-t[i])) * (1 + A)
      B = exp(-delta*(t[i+1]-t[i])) * B + sum(exp(-delta*B_sequence[[i+1]]))
    }
  }
  return(-logl)
}

## Function to calculate B_sequence
B_sequence_calculator = function(t, t_prime){
  B_sequence = list()
  index = ifelse(t_prime[1] < t[1], max(which(t_prime < t[1])), 0)
  if(index != 0){
    B_sequence[[1]] = t[1] - t_prime[1:index]
  } else {
    B_sequence[[1]] = numeric(0)
  }
  for(i in 2:length(t)){
    B_sequence[[i]] = t[i] - t_prime[t_prime < t[i] & t_prime >= t[i-1]]
  }
  return(B_sequence)
}

## Optimise the log-likelihood
optim_logl = function(starting_values,t,t_prime,T,B_sequence){
  neg_logl = function(x) log_likelihood(params=x, t=t, t_prime=t_prime, T=T, B_sequence=B_sequence)
  return(optim(par=starting_values, fn=neg_logl))
}

## Load some utility functions
source('hawkes_double.R')
## If in hawkes_double.R you change test = FALSE to test = TRUE on line 76, you will get some 
## simulated data with parameters specified in the file hawkes_double.R 
## The event times are in t and t_prime ** which would replace the Hyde Park event times! **

## Calculate B_sequence
params = c(0.1, 0.25, 1.0, 0.1, 0.25)
B_seq = B_sequence_calculator(t, t_prime)

## MLE
res = optim_logl(starting_values=transform_parameters(params, to_constrained=FALSE),
                 t=t, t_prime=t_prime, T=T, B_sequence=B_seq)
p = transform_parameters(res$par, to_constrained=TRUE)

## Constrained optimisation
fixed = c(-Inf,0)
optim_logl_constrained = function(starting_values,t,t_prime,T,B_sequence){
  neg_logl = function(x) log_likelihood(params=c(x,fixed), t=t, 
                                        t_prime=t_prime, T=T, B_sequence=B_sequence)
  return(optim(par=starting_values, fn=neg_logl))
}
res_model1 = optim_logl_constrained(starting_values=transform_parameters(params, to_constrained=FALSE)[1:3],
                 t=t, t_prime=t_prime, T=T, B_sequence=B_seq)
p2 = transform_parameters(c(res_model1$par,fixed), to_constrained=TRUE)

## Important to check the dependence on starting values --> better results here!!!
res = optim_logl(starting_values=rep(0,5),
                 t=t, t_prime=t_prime, T=T, B_sequence=B_seq)
p = transform_parameters(res$par, to_constrained=TRUE)
## If simulated data are used, something very close to the 'truth' is obtained here 
## params = c(0.1, 0.25, 1.0, 0.1, 0.25) was used in the simulation 

res_model1 = optim_logl_constrained(starting_values=c(0,0,0),
                                    t=t, t_prime=t_prime, T=T, B_sequence=B_seq)
p2 = transform_parameters(c(res_model1$par,fixed), to_constrained=TRUE)

## Other model: only consider the arrival times of the journeys ending at station i
optim_logl_constrained2 = function(starting_values,t,t_prime,T,B_sequence){
  neg_logl = function(x) log_likelihood(params=c(x[1],fixed,x[2],x[3]), t=t, 
                                        t_prime=t_prime, T=T, B_sequence=B_sequence)
  return(optim(par=starting_values, fn=neg_logl))
}
res_model2 = optim_logl_constrained2(starting_values=c(0,0,0),
                                    t=t, t_prime=t_prime, T=T, B_sequence=B_seq)
p3 = transform_parameters(c(res_model2$par,fixed), to_constrained=TRUE)[c(1,4,5,2,3)]
