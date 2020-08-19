source('likelihood_func.R')

theta <- c(1,0.00166,5,0.5,0.05)

negative_ll <- function(theta,inddata,outdata){
  -LogLikelihood(theta,inddata,outdata)
}

result <- optim(par = theta,fn = negative_ll, inddata = inddata, outdata = outdata, control=list(trace=1))

# Plot the optimisation result
res_multi_seird <- ForwardSimulation(result$par)
plot_results(res_multi_seird, inddata, outdata)

