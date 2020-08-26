# mcmc mtrop
library('mcmc')
library('readr')
library(deSolve) 
library(ggplot2)
library(lubridate)
library(jsonlite)
library(tidyverse)
library(dplyr)
library(grid)
library(readr)
library(scales)

ddata <- read_excel("Desktop/Cov model/uk daily.xlsx")   # the outside death population
data <- read_excel("Desktop/Cov model/0306-0612care_home_deaths.xlsx")  # the death population inside care homes
outdata <- ddata$`deaths outside the care homes`  # outside the carehomes
inddata <- data$`Deaths involving COVID-19 of care home residents` # inside
LogLikelihood <- function(inddata, outdata)function(theta){
  penalty_for_negative_params <- 0
  if(theta[1]<0)
  {
    penalty_for_negative_params <- penalty_for_negative_params -1e9*(1-theta[1])
  }
  if(theta[2]<0)
  {
    penalty_for_negative_params <- penalty_for_negative_params -1e9*(1-theta[2])
  }
  if(theta[3]<0)
  {
    penalty_for_negative_params <- penalty_for_negative_params -1e9*(1-theta[3])
  }
  if(theta[4]<0)
  {
    penalty_for_negative_params <- penalty_for_negative_params -1e9*(1-theta[4])
  }
  if(theta[5]<0)
  {
    penalty_for_negative_params <- penalty_for_negative_params -1e9*(1-theta[5])
  }
  
  if (penalty_for_negative_params < 0){
    return(penalty_for_negative_params)
  }
  
  res_multi_seird <- ForwardSimulation(theta)
  
  simulate_death <- res_multi_seird$Dc
  outside_death <- res_multi_seird$Dn
  daily_inside_death <- simulate_death[2:100]-simulate_death[1:99]
  daily_outside_death <- outside_death[2:100]-outside_death[1:99]
  
  iout <- ll(inddata,daily_inside_death) # the loglikelihood values for inside death
  oout <- ll(outdata,daily_outside_death) # the loglikelihood values for outside death
  return(oout + iout)
}

# loglikelihood function
ll <- function(x,mu){
  sigma <- 100
  #theta is a vector containing the two parameters of interest
  n <- length(x)
  logl <- -.5*n*log(2*pi) -.5*n*log(sigma^2)-
    (1/(2*(sigma^2)))*(sum((x - mu)^2))
  return(logl)
}

ForwardSimulation <- function(theta){
  Nn <- 64582000 # total population outside the care homes 
  Nc <- 418000 # total population inside the care homes of England and Wales
  In0 <- theta[3] # initial number of infectious population outside the care homes
  Ic0 <- 0 # initial number of infectious population inside care homes
  En0 <- 0 # initial number of exposed population
  Rn0 <- 0 # initial number of recovered population
  Dn0 <- 0 # initial number of death population
  Ec0 <- 0 # initial number of exposed population
  Rc0 <- 0 # initial number of recovered population
  Dc0 <- 0 # initial number of death population
  Sn0 = (Nn - In0 - Rn0 -  Dn0) # initial number of susceptible population
  Sc0 = (Nc - Ic0 - Rc0 -  Dc0) # initial number of susceptible population inside care homes
  init <- c(Sn = Sn0, En = En0, In = In0, Rn = Rn0, Dn = Dn0, Nn = Nn,
            Sc = Sc0, Ec = Ec0, Ic = Ic0, Rc = Rc0, Dc = Dc0, Nc = Nc)	
  time <- seq(0, 99, 1) 
  
  pars<-c( 
    betann = theta[1],	# infection rate for Sn to In
    betanc = theta[1],	# infection rate for Sn to Ic
    betacn = theta[1],	# infection rate for Sc to In
    betacc = theta[2],	# infection rate for Sc to Ic
    alphan = 0.8 ,	# transmission rate from exposed individual to infected individual
    gamman = 0.65,	# recovery rate 
    mun = theta[5],	# death rate outside the care homes
    alphac = 0.8,	# transmission rate from exposed individual to infected individual inside care homes
    gammac = 0.65,	# recovery rate inside care homes
    muc = theta[4]	# death rate inside care homes
  ) 
  res_multi_seird <- as.data.frame(lsoda(y = init, times = time, func = multi_seird, parms = pars))
  return(res_multi_seird)
}

# define the cofficient g(t) for death rate 
# Define the parameters
go <- 34.84 # 40
gb <- 0.21 # lockdown_baseline
gr <- 0.15 # rate_of_lockdown
go <- 34.84 # lockdown_offsite
gf <- 0.00606 #lockdown_fatigue_rate
g <- function(p,i){
  if (i >= go)
  {fi <- -(1 - gb)*exp(- (i - go)*gf)+(1 - gb)}
  else
  {fi <- 0}
  gi <- 1 - p*(0.5*(1 - gb)*(tanh(gr*(i - go))+1 )+ fi)
  return(gi)
}
multi_seird <- function(time, state, pars){ 
  with(as.list(c(state, pars)),{ 
    dSn <- -Sn * (betann*g(1,time)*In + betanc*g(0.5,time)*Ic)/Nn
    dEn <- Sn * (betann*g(1,time)*In + betanc*g(0.5,time)*Ic)/Nn - En * alphan 
    dIn <- En * alphan -  In * gamman
    dRn <- In * (gamman - mun)
    dDn <- In * mun
    dNn <- dSn + dEn +dIn + dRn + dDn 
    dSc <- -Sc * (betacn*g(0.5,time)*In + betacc*Ic)/Nc
    dEc <- Sc * (betacn*g(0.5,time)*In + betacc*Ic)/Nc - Ec * alphac 
    dIc <- Ec * alphac - Ic * gammac
    dRc <- Ic * (gammac - muc)
    dDc <- Ic * muc
    dNc <- dSc + dEc + dIc + dRc + dDc
    return(list(c(dSn,dEn,dIn,dRn,dDn,dNn,
                  dSc,dEc,dIc,dRc,dDc,dNc)))
  })
} 
# begin mcmc
set.seed(46) # to get reproducible results
theta.init <- c(1,0.00166,25,0.5,0.05)
model <- LogLikelihood(inddata = inddata, outdata = outdata)
out <- metrop(model, theta.init, 1e3)
names(out)
out$accept
out <- metrop(out, scale = 0.001)
out$accept
# [1] 0.676
out <- metrop(out, scale = 0.005)
out$accept
# [1] 0.203

# 3 disanostics
out <- metrop(out, nbatch = 1e4)
t.test(out$accept.batch)$conf.int
out$time

plot(ts(out$batch))
acf(out$batch)

# 4 Monte Carlo Estimates and Standard Errors
out <- metrop(out, nbatch = 100, blen = 100,
            outfun = function(z) c(z, z^2))
t.test(out$accept.batch)$conf.int
out$time

apply(out$batch, 2, mean)
foo <- apply(out$batch, 2, mean)
mu <- foo[1:5]
sigmasq <- foo[6:10] - mu^2
mu
sigmasq

mu.mcse <- apply(out$batch[ , 1:5], 2, sd) / sqrt(out$nbatch)
mu.mcse

#4.2 Functions of Means
u <- out$batch[ , 1:5]
v <- out$batch[ , 6:10]
ubar <- apply(u, 2, mean)
vbar <- apply(v, 2, mean)
deltau <- sweep(u, 2, ubar)
deltav <- sweep(v, 2, vbar)
foo <- sweep(deltau, 2, ubar, "*")
sigmasq.mcse <- sqrt(apply((deltav - 2 * foo)^2, 2, mean) / out$nbatch)
sigmasq.mcse
sqrt(mean(((v[ , 2] - vbar[2]) - 2 * ubar[2] * (u[ , 2] - ubar[2]))^2) /
       out$nbatch)
#4.3 Functions of Functions of Means
sigma <- sqrt(sigmasq)
sigma.mcse <- sigmasq.mcse / (2 * sigma)
sigma
sigma.mcse

#5 A Final Run
out <- metrop(out, nbatch = 500, blen = 400)
t.test(out$accept.batch)$conf.int
out$time
foo <- apply(out$batch, 2, mean)
mu <- foo[1:6]
sigmasq <- foo[6:10] - mu^2
mu
sigmasq
mu.mcse <- apply(out$batch[ , 1:5], 2, sd) / sqrt(out$nbatch)
mu.mcse
u <- out$batch[ , 1:5]
v <- out$batch[ , 6:10]
ubar <- apply(u, 2, mean)
vbar <- apply(v, 2, mean)
deltau <- sweep(u, 2, ubar)
deltav <- sweep(v, 2, vbar)
foo <- sweep(deltau, 2, ubar, "*")
sigmasq.mcse <- sqrt(apply((deltav - 2 * foo)^2, 2, mean) / out$nbatch)
sigmasq.mcse
sigma <- sqrt(sigmasq)
sigma.mcse <- sigmasq.mcse / (2 * sigma)
sigma
sigma.mcse


