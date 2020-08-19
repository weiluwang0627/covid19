# plot of the death population
# import the dataset
ddata <- read_excel("Desktop/uk daily.xlsx")   # the total death population
data <- read_excel("Desktop/0306-0612care_home_deaths.xlsx")  # the death population inside care homes
dailydata <- ddata$`Daily change in deaths`
inddata <- data$`Deaths involving COVID-19 of care home residents`
date <- data$Date
datet <- ddata$Date
timei <- seq(1,99,1)
inside <- data.frame(timei,date,inddata)
totald <- data.frame(timei,date,ddata)
model <-data.frame(date,inside,totald)
ggplot(model,aes(date)) + 
  geom_line(aes(x = date, y = dailydata, col = 'total daily death population'))+
  geom_line(aes(x = date, y = inddata, col='daily death population inside care homes'))+
  scale_colour_manual("",
                      values=c(
                        "daily death population inside care homes" = "forestgreen", "total daily death population" = "cornflowerblue" 
                      ) ) + 
  scale_y_continuous('') + labs(title = "Death Population") + scale_x_date(breaks=as.date(c("2020-03-06","2020-04-06","2020-05-06"))


# define the normal loglikelihood function for total death
tnll <- function(x){
  mu1 = 100  # mu1 here is an estimate parameter 
  sigma = 100
  #theta is a vector containing the two parameters of interest
  n <- seq(1,146,1)
  logl <- -.5*n*log(2*pi) -.5*n*log(sigma^2)-
    ((1/(2*(sigma^2)))*(sum((x - mu1)^2))
     return(-logl)
} 

tout <- tnll(dailydata) # the loglikelihood values for total death

inll <- function(x){
  mu1 = 100
  sigma = 100
  #theta is a vector containing the two parameters of interest
  n <- seq(1,167,1)
  logl <- -.5*n*log(2*pi) -.5*n*log(sigma^2)-
    ((1/(2*(sigma^2)))*(sum((x - mu1)^2))
     return(-logl)
} 
iout <- inll(inddata)
out <- iout + tout
out  # the sum of them
# find the MLE
max(out)