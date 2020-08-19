library(deSolve) 
library(ggplot2)
library(lubridate)
library(jsonlite)
library(tidyverse)
library(dplyr)
library(grid)
library(readr)
library(scales)


# import the ONS dataset of deaths
data <- read_excel("Desktop/death_data.xlsx")
view(data)

cdeath <- data$`ONS Deaths involving COVID-19`# death population caused by covid_19
total <- data$`ONS All deaths year-to-date` # total deaths of this year
ldeath <- data$`ONS 2019` # death population of 2019 
time <- as.Date(data$Date)
days <- seq(1,168,1)
# plot 
df = data.frame(days,ldeath,cdeath,total)

# build the SEIRD model to get the simulate death population number
seird <- function(time, state, pars){ 
  with(as.list(c(state, pars)),{ 
    dS <- -S * beta * I/N 
    dE <- S * beta * I/N - E * alpha 
    dI <- E * alpha - I * gamma
    dR <- I * (gamma - mu)
    dD <- I * mu
    dN <- dS + dE + dI + dR + dD 
    return(list(c(dS,dE,dI,dR,dD,dN)))
  })
} 
# define the parameters inside care homes
N <- 418000 # total population inside care homes
I0 <- 2 # initial number of infected population
E0 <- 0 # initial number of exposed population
R0 <- 0 # initial number of recovered population
D0 <- 0 # initial number of death population
S0 = N - I0 - R0 -  D0 # initial number of susceptible population
init <- c(S = S0, E = E0, I = I0, R = R0, D = D0 , N = N)	
time <- seq(1, 168, 1)  
pars<-c( 
  beta = 0.376,	# infection rate 0.386
  alpha = 1,	# transmission rate from exposed individual to infected individual
  gamma = 0.25,	# recovery rate 
  mu= 0.0202	# death rate
) 
y1 <- data.frame(time,cdeath)
y2 <- as.data.frame(lsoda(y = init, times = time, func = seird, parms = pars)) 
df <- data.frame(time,y1,y2)
# plot the comparison of the real death population and simulation of death population
ggplot(df,aes(time)) + 
  geom_line(aes(x = time, y = D, col = 'Simulate Death Population'))+
  geom_line(aes(x = time, y = cdeath, col='Death Population Caused by COVID-19'))+
scale_colour_manual("",
                    values=c(
                      "Simulate Death Population" = "forestgreen", "Death Population Caused by COVID-19" = "cornflowerblue" 
                    ) 
) + scale_y_continuous('')+ labs(title = "Death Population Comparisions")

# find the cofficient parameter g(t) for death rate mu
# first values of the death population are zeros
# we take the inverse of them to calculate the rate
g <- ((y1$cdeath)/(y2$D))
gt <- data.frame(time,g)
ggplot(gt,aes(time))+
  geom_line(aes(x = time, y = g, col = 'Cofficient of The Death Rate'))+
  scale_colour_manual("",
                      values=c(
                        "Cofficient of The Death Rate" = "red"))+
                        scale_y_continuous('')+ 
                        labs(title = "The Cofficient Parameter g(t) for Death Rate")



