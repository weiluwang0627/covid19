# This file is just an example call with a given set of parameters which
# we will use to specify a likelihood value
source('likelihood_func.R')

# Define the parameters
theta <- c(1,1,2,0.064,0.02)

res_multi_seird <- ForwardSimulation(theta)

plot_results(res_multi_seird, inddata, outdata)



