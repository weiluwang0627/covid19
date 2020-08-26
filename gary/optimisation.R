setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('likelihood_func.R')

#theta <- c(4000, 1.1, 232)
negative_ll_just_outside <- function(theta,inddata,outdata){
  -LogLikelihoodJustOutside(theta,inddata,outdata)
}

theta <- c(3800, 1.1, 0.01, 0.01 , 1.1, 0.04, 200, 35)
negative_ll <- function(theta,inddata,outdata){
  -LogLikelihood(theta,inddata,outdata)
}

result <- optim(par = theta,fn = negative_ll, inddata = inddata, outdata = outdata, control=list(trace=1))

# Plot the optimisation result
par(mfrow=c(3,1)) 
res_multi_seird <- ForwardSimulation(result$par)
plot_results(res_multi_seird, inddata, outdata)

