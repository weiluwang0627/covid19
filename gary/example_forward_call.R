# This file is just an example call with a given set of parameters which
# we will use to specify a likelihood value
source('likelihood_func.R')

# Define the parameters
theta <- c(4170, 1.1, 0.0026, 0.045 , 1.4, 0.026, 197, 45)

res_multi_seird <- ForwardSimulation(theta)

par(mfrow=c(3,1)) 
plot_results(res_multi_seird, inddata, outdata)



