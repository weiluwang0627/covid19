setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('likelihood_func.R')

#install.packages('adaptMCMC')
require('adaptMCMC')

# Result of an optimization run
theta <- c(4170, 1.1, 0.0026, 0.045 , 1.4, 0.026, 197, 45)

num_samples <- 200000
result <- MCMC(p = LogLikelihood, init = theta, scale = 0.01*theta, n=num_samples, adapt=TRUE, acc.rate=0.2, inddata = inddata, outdata = outdata)

# Plot the MCMC results
result$acceptance.rate

plot(ts(result$samples))

par(mfrow=c(3,2)) 
hist(result$samples[,1])
hist(result$samples[,2])
hist(result$samples[,3])
hist(result$samples[,4])
hist(result$samples[,5])

dev.new()
max_posterior = which.max(result$log.p)
sim_max_posterior <- ForwardSimulation(result$samples[max_posterior,])
par(mfrow=c(3,1)) 
plot_results(sim_max_posterior, inddata, outdata)

for (i in seq(from = num_samples/2, to = num_samples, by = 500))
{
  sim_sample <- ForwardSimulation(result$samples[i,])
  plot_results(sim_sample, inddata, outdata)
}
