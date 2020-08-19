# mcmc
library('mcmc')
# Define the parameters
go <- 34.84
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
# Define the parameters
# unif parameters
theta.init <- c(1*exp(-5),1*exp(-5),1*exp(-5),1*exp(-3),1*exp(-3),64582000,418000)
theta[1] <- runif(1,1*exp(-5),1)
theta[2] <- runif(1,1*exp(-5),1)
theta[3] <- runif(1,1*exp(-5),1)
theta[4] <- runif(1,1*exp(-3),1)
theta[5] <- runif(1,1*exp(-3),1)
theta[6] <- runif(1,64582000,64582100)
theta[7] <- runif(1,418000,418100)
Nn <- 64582000 # total population outside the care homes 
Nc <- 418000 # total population inside the care homes of England and Wales
In0 <- theta[3] # initial number of infectious population outside the care homes
Ic0 <- theta[3] # initial number of infectious population inside care homes
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
time <- seq(1, 99, 1) 
i <- time
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

# import the dataset
ddata <- read_excel("Desktop/uk daily.xlsx")   # the outside death population
data <- read_excel("Desktop/0306-0612care_home_deaths.xlsx")  # the death population inside care homes
outdata <- ddata$`deaths outside the care homes`  # outside the carehomes
inddata <- data$`Deaths involving COVID-19 of care home residents` # inside

model <- function(x,mu)function(theta){
  for (i in 2:99){ 
  death_llh <- as.data.frame(lsoda(y = init,
            times = i, func = multi_seird, parms = pars))
  sigma <- 100
  #theta is a vector containing the two parameters of interest
  n <- seq(1,99,1)
  logl <- -.5*i*log(2*pi) -.5*i*log(sigma^2)-
    (1/(2*(sigma^2)))*(sum((x - mu)^2))
  return(-logl)}
}
inside_simulate_death <- (as.data.frame(death_llh))$Dc
outside_simulate_death <- (as.data.frame(death_llh))$Dn
mu3 <- inside_simulate_death[2:100]-inside_simulate_death[1:99]
mu4 <-  outside_simulate_death[2:100]-outside_simulate_death[1:99]

inside_out <- model(inddata,mu3)
outside_out <- model(outdata,mu4)
total_out <- inside_out + outside
#2 beginning mcmc
set.seed(42) # to get reproducible results
theta.init <- c(1*exp(-5),1*exp(-5),1*exp(-5),1*exp(-3),1*exp(-3))
out <- metrop(inside_out,theta.init, 1e3)
names(out)




