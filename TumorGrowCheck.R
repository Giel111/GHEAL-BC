library(tidyverse)
library(ggplot2)
library(gridExtra)

hour = 1;
minute = hour/60;
second=minute/60;
day=hour*24;
year = day*365; 
month = year/12 #not exactly but fine


max_size = 80
start_size = 0.25
Vc = (4/3)*pi*(start_size/2)**3 #Volume at start
Vm = (4/3)*pi*(max_size/2)**3 #Max volume
S = 5*10e-11
G = 1.07/month #0.341/62 # per day
sigma = 1.31 # (0.539-0.133)/3.92 * np.sqrt(43) 


GompGrow <- function(age,onsetage,GR,IsCured){
  if(IsCured ==1){return(0)} else{
    t <- (age-onsetage ) / month
    if (t<0){
      size<-0
    } else {
      #Determine size of tumour MANC
      Volume <- Vm/(1+((Vm/Vc)^0.25-1)*exp(-0.25*GR*t))^4 #tumour volume at time t
      size <- 2*(Volume/(4/3*pi))^(1/3)
    }
    return(size)
  }
}

FindTimeAtSize <- function(size, onsetage, GR, max_size) {
  
  # Solve for the time delta using a binary search algorithm
  low <- 0
  high <- 1000  # set an arbitrarily high upper bound for time delta
  while (high - low > 1e-6) {
    mid <- (low + high) / 2
    #size_mid <- GompGrow(onsetage + mid, onsetage, GR, IsCured=0)
    Volume_mid <- Vm/(1+((Vm/Vc)^0.25-1)*exp(-0.25*GR*mid))^4 #tumour volume at time t
    size_mid <- 2*(Volume_mid/(4/3*pi))^(1/3)
    if (size_mid < size) {
      low <- mid
    } else {
      high <- mid
    }
  }
  
  # Return the estimated time when the given size was reached
  return(onsetage + high)
}

GompGrow2 <- function(age,onsetage,GR,IsCured,max_size,RegressionSize,StagnateSize){
  Vm = (4/3)*pi*(max_size/2)**3 
  size <- 0
  delta = age-onsetage
  if (delta<0){
    size <-0
    return (size)
  } 
  #Determine size of tumour MANC
  Volume <- Vm/(1+((Vm/Vc)^0.25-1)*exp(-0.25*GR*delta))^4 #tumour volume at time t
  size <- 2*(Volume/(4/3*pi))^(1/3)
  if (size>StagnateSize){
    size <- StagnateSize
  }
  if (size>RegressionSize){
    size <- RegressionSize
    time_reached = FindTimeAtSize(RegressionSize,onsetage,GR,max_size)
    delta2 = delta-time_reached
    size2 = size - (size*pnorm(delta2,50,25))
    return(size2)
  }
  return(size)
}


#This function uses the same GompGrow2 function as in the TimeDelta function to estimate the tumor size at each time point, and performs a binary search to

# Create a sequence of numbers between -10 and 10 incrementing by 0.2.
x <- seq(0,100,by = .2)
# Choose the mean as 2.5 and standard deviation as 2. 
y <- -1*pnorm(x, mean = 50, sd = 25)
plot(x,y)


# Create a vector of time points from 0 to 1000
time_points <- seq(95,95 + 12*6, by = 0.1)

Grow_Gamma_mean <- 0.12
Grow_Gamma_sd <- 0.012

# Create a vector of the corresponding cell sizes at each time point
GR1 = 0.25# rlnorm(1,meanlog=Grow_log_norm_mean,sdlog=Grow_log_norm_sd)
GR1 = rgamma(1,Grow_Gamma_mean,Grow_Gamma_sd)
GR1 = qgamma(0.95,Grow_Gamma_mean,Grow_Gamma_sd)
GR1 = qgamma(0.05,Grow_Gamma_mean,Grow_Gamma_sd)
# Calculate the median growth rate that corresponds to 84 months
GR50 <- log(2)/(84/12)

# Calculate the standard deviation that corresponds to the 95th and 25th percentile values
q95 <- qnorm(0.95, log(GR50), log(144/24))
q25 <- qnorm(0.25, log(GR50), log(144/24))
sd_log <- (q95 - q25) / (2 * qnorm(0.75))

# Generate 1000 random samples from the lognormal distribution
#set.seed(123) # for reproducibility
GR1 <- rlnorm(n = 1, meanlog = log(GR50), sdlog = sd_log)

# Set the percentiles and target time values
p1 <- 0.95
p2 <- 0.25
t1 <- 6
t2 <- 48

# Calculate the shape and rate parameters of the gamma distribution
q1 <- qgamma(p1, shape = 2, rate = t1/2)
q2 <- qgamma(p2, shape = 2, rate = t2/2)

# Calculate the mean and standard deviation of the distribution
mean_gr <- (q1 + q2) / 2
sd_gr <- (q2 - q1) / (2 * qnorm(p1 + p2 - 1) * sqrt(2))

# Generate GR values from the gamma distribution
#set.seed(123)  # for reproducibility
GR1 <- rgamma(1, shape = (mean_gr / sd_gr)^2, rate = mean_gr / sd_gr^2)

print(GR1)
cell_sizes <- sapply(time_points, function(t) GompGrow2(t,onsetage=0,GR=GR1,IsCured=0,max_size=128,RegressionSize = 1000,StagnateSize = 1000))

# Plot the cell sizes over time
plot(time_points, cell_sizes, type = "l", xlab = "Time (months)", ylab = "Tumor diameter (mm)",ylim=c(0,150))


vect <- rlnorm(10000,meanlog=Grow_log_norm_mean,sdlog=Grow_log_norm_sd)
vect<- sapply(vect,function(a) min(25,a))
vect <- rgamma(1000,Grow_Gamma_mean,Grow_Gamma_sd) 
vect<- sapply(vect,function(a) min(25,a))
hist(vect)

library(tidyverse)
# Create a vector of time points from 0 to 1000
time_points <- seq(-5, 0 + 12 * 13, by = 0.1)

# Generate gamma distribution parameters based on mean and standard deviation of growth rates
shape <- (mean_gr / sd_gr)^2
rate <- mean_gr / sd_gr^2

# Generate growth rates for each quantile
growth_rates <- qgamma(c(0.05, 0.25, 0.5, 0.75, 0.95), shape = shape, rate = rate)

# Calculate cell sizes for each growth rate at each time point
cell_sizes <- sapply(growth_rates, function(gr) {
  sapply(time_points, function(t) {
    GompGrow2(t, onsetage = 0, GR = gr, IsCured = 0, max_size = 128,
              RegressionSize = 1000, StagnateSize = 1000)
  })
})

# Plot the cell sizes over time for each growth rate
plot(time_points, cell_sizes[,1 ], type = "l", xlab = "Time (months)", ylab = "Cell size (mm)", col = "red",ylim=c(0,150))
lines(time_points, cell_sizes[,2 ], col = "orange")
lines(time_points, cell_sizes[,3 ], col = "blue")
lines(time_points, cell_sizes[,4 ], col = "green")
lines(time_points, cell_sizes[,5 ], col = "purple")

# Add legend
legend(x=100,y=100, legend = c("95th percentile", "75th percentile", "Median", "25th percentile",
                              "5th percentile"), col = c("purple", "green", "blue", "orange", "red"),
       lty = 1,cex=0.75)

## CONCLUDING:

# now we've fond the right distributino and the right params:
# gamma dist
# shape = 1.567742
# rate = 1.933883
