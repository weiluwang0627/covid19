setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('likelihood_func.R')

#install.packages('adaptMCMC')
require('adaptMCMC')

# Result of an optimization run
theta <- c(0.5022554, 0.5236479, 20.3501919,  0.2732997,  0.4517723)

result <- MCMC(p = LogLikelihood, init = theta, scale = 0.01*theta, n=50000, adapt=TRUE, acc.rate=0.234, inddata = inddata, outdata = outdata)

# Plot the MCMC results
result$acceptance.rate

plot(ts(result$samples))

par(mfrow=c(3,2)) 
hist(result$samples[,1])
hist(result$samples[,2])
hist(result$samples[,3])
hist(result$samples[,4])
hist(result$samples[,5])
