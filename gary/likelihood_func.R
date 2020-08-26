# Create a ForwardSimulation method 
# build the SEIRD model
# define the cofficient g(t) for death rate 
# Define the parameters

#install.packages('deSolve')
#install.packages('readxl')

require('deSolve')
require('readxl')

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# import the dataset (outside functions or loops)
data <- read_excel("../0306-0612care_home_deaths.xlsx")  # the death population inside care homes
inddata <- data$`Deaths involving COVID-19 of care home residents` # inside

ddata <- read_excel("../uk daily.xlsx")   # the outside death population
outdata <- ddata$`deaths outside the care homes`  # outside the carehomes

stopifnot(length(inddata)==length(outdata))
days_to_simulate <- length(inddata)
num_days_infectious_to_death <- 15

LogLikelihoodJustOutside <- function(theta, inddata, outdata){
  theta_larger <- c(theta[1:2],0,0,0,0,theta[3],10)
  LogLikelihood(theta=theta_larger,inddata,outdata)
}

LogLikelihood <- function(theta, inddata, outdata){
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
  daily_inside_death <- simulate_death[2:(days_to_simulate+1)]-simulate_death[1:days_to_simulate]
  daily_outside_death <- outside_death[2:(days_to_simulate+1)]-outside_death[1:days_to_simulate]
  
  # Now account for the delay from leaving infectious group to death occurring
  daily_inside_death <- c(numeric(num_days_infectious_to_death),daily_inside_death[1:(days_to_simulate-num_days_infectious_to_death)])
  daily_outside_death <- c(numeric(num_days_infectious_to_death),daily_outside_death[1:(days_to_simulate-num_days_infectious_to_death)])
  
  iout <- ll(inddata,daily_inside_death, theta[7]) # the loglikelihood values for inside death
  oout <- ll(outdata,daily_outside_death, theta[8]) # the loglikelihood values for outside death
  return(oout + iout)
}

# loglikelihood function
ll <- function(x,mu,sigma){
  #theta is a vector containing the two parameters of interest
  n <- length(x)
  logl <- -.5*n*log(2*pi) -.5*n*log(sigma^2)-
    (1/(2*(sigma^2)))*(sum((x - mu)^2))
  return(logl)
}

ForwardSimulation <- function(theta){
  Nn <- 64582000 # total population outside the care homes 
  Nc <- 418000 # total population inside the care homes of England and Wales
  In0 <- theta[1] # initial number of infectious population outside the care homes
  Ic0 <- 0 # initial number of infectious population inside care homes
  En0 <- 0 # initial number of exposed population
  Rn0 <- 0 # initial number of recovered population
  Dn0 <- 0 # initial number of death population
  Ec0 <- 0 # initial number of exposed population
  Rc0 <- 0 # initial number of recovered population
  Dc0 <- 0 # initial number of death population
  Sn0 = (Nn - In0 - Rn0 -  Dn0) # initial number of susceptible population
  Sc0 = (Nc - Ic0 - Rc0 -  Dc0) # initial number of susceptible population inside care homes
  init <- c(Sn = Sn0, En = En0, In = In0, Rn = Rn0, Dn = Dn0,
            Sc = Sc0, Ec = Ec0, Ic = Ic0, Rc = Rc0, Dc = Dc0)	
  time <- seq(0, days_to_simulate, 1) 
  
  pars<-c( 
    betann = theta[2],	# infection rate for Sn to In
    betanc = theta[3],	# infection rate for Sn to Ic
    betacn = theta[4],	# infection rate for Sc to In
    betacc = theta[5],	# infection rate for Sc to Ic
    alphan = 0.2 ,  # transmission rate from exposed individual to infected individual
    gamman = 0.5,	  # removal rate 
    mun = 0.5/100,	# death rate outside the care homes ~ 1% IFR
    alphac = 0.2,	  # transmission rate from exposed individual to infected individual inside care homes
    gammac = 0.5,	  # removal rate inside care homes
    muc = theta[6],	# death rate inside care homes ~ needs to be much higher than 1% IFR
    Nn = Nn,
    Nc = Nc
  ) 
  res_multi_seird <- as.data.frame(lsoda(y = init, times = time, func = multi_seird, parms = pars))
  return(res_multi_seird)
}

go <- 34.84
gb <- 0.21 # lockdown_baseline
gr <- 0.15 # rate_of_lockdown
go <- 34.84 # lockdown_offset
gf <- 0.00606 #lockdown_fatigue_rate
g <- function(p,i){
  if (i >= go)
  {fi <- -(1 - gb)*exp(- (i - go)*gf)+(1 - gb)}
  else
  {fi <- 0}
  gi <- 1 - p*(0.5*(1 - gb)*(tanh(gr*(i - go))+1 )- fi)
  return(gi)
}

multi_seird <- function(time, state, pars){ 
  with(as.list(c(state, pars)),{ 
    dSn <- -Sn * (betann*g(1,time)*In + betacn*g(0.5,time)*Ic)/Nn
    dEn <- Sn * (betann*g(1,time)*In + betacn*g(0.5,time)*Ic)/Nn - En * alphan 
    dIn <- En * alphan -  In * gamman
    dRn <- In * (gamman - mun)
    dDn <- In * mun
    dSc <- -Sc * (betanc*g(0.5,time)*In + betacc*g(0.5,time)*Ic)/Nc
    dEc <- Sc * (betanc*g(0.5,time)*In + betacc*g(0.5,time)*Ic)/Nc - Ec * alphac 
    dIc <- Ec * alphac - Ic * gammac
    dRc <- Ic * (gammac - muc)
    dDc <- Ic * muc
    return(list(c(dSn,dEn,dIn,dRn,dDn,
                  dSc,dEc,dIc,dRc,dDc)))
  })
} 

plot_results = function(res_multi_seird, inddata, outdata){
  simulate_death <- res_multi_seird$Dc
  outside_death <- res_multi_seird$Dn
  daily_inside_death <- simulate_death[2:(days_to_simulate+1)]-simulate_death[1:days_to_simulate]
  daily_outside_death <- outside_death[2:(days_to_simulate+1)]-outside_death[1:days_to_simulate]
  
  # Now account for the delay from leaving infectious group to death occurring
  daily_inside_death <- c(numeric(num_days_infectious_to_death),daily_inside_death[1:(days_to_simulate-num_days_infectious_to_death)])
  daily_outside_death <- c(numeric(num_days_infectious_to_death),daily_outside_death[1:(days_to_simulate-num_days_infectious_to_death)])
  
  days <- 1:days_to_simulate
  
  par(mfg=c(1,1))
  plot(days,inddata)
  lines(days,daily_inside_death)
  
  par(mfg=c(2,1))
  plot(days,outdata)
  lines(days,daily_outside_death)
  
  par(mfg=c(3,1))
  google <- vector()
  for (i in days)
  {
      google <- append(google,g(1,i))
  }
  plot(days, google)
}
