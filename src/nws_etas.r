################################################################################
# 
# Fit the ETAS function to the earthquake occurrence data, using the  
# Saichev & Sornette (2006) approach
#
# Stefan Revets
#
# 2012-03-28 Version 1.0
#
################################################################################

# Define the ETAS function
ETAS <- function(x, alpha, n, theta, rho){
  (alpha * n * theta * rho^theta * x^(-1-theta) + (1 - n + alpha * n * rho^theta * x^(-theta))^2) * exp(-(1 - n) * x - alpha * n *rho^theta * x^(1-theta) / (1-theta))
  }

Corral <- function(x, C, alpha, g){
  C * exp(-(x/alpha)) * (x/alpha)^(g-1) / (alpha * gamma(g))
  }

# Read the data from file. Expected format is:
# UTC(YYYY-MM-DD hh:mm:ss[.ss])  Magnitude Latitude  Longitude Depth
# preferably in ascending UTC order
all_events <- read.table("nws_data.txt", header = T, sep = "\t")
# and convert the UTC time into (decimal) days
all_events$UTC <- as.numeric(difftime(all_events$UTC,all_events$UTC[1], units="days"))
# There is a problem with the magnitudes pre-1992 (see Sagar & Leonard) for
# NW Australia: these are about 0.5 too high
all_events$Magnitude[1:364] <- all_events$Magnitude[1:364] - 0.5

## Testing sensitivity to magnitude estimate errors
## Just perturb the actual magnitudes by e(0.0, 0.2)
#My_perturb <- rnorm(length(all_events$Magnitude), mean=0.0, sd = 0.2)
#all_events$Magnitude <- all_events$Magnitude + My_perturb
all_events <- all_events[3:length(all_events$UTC),]
#all_events <- all_events[95:580,]
Max_Timespan <- round(max(all_events$UTC) - min(all_events$UTC))

# Set up the time bins used to calculate the event densities
My_Breaks <- c(0,Max_Timespan*2^-seq(15,0))

# Set up the holder for the processed data
Waiting_Data <- data.frame(Time=rep(0,times=112) ,Density=rep(0,times=112))

# Set lower magnitude cut-off, and  readjust the database accordingly
for (Counter in 1:7)
  {
  Min_Mag <-  1.5 + Counter * 0.5
  events <- subset(all_events, Magnitude > Min_Mag, select = c(UTC,Magnitude,Longitude,Latitude,Depth))
  
# Get the average seismic rate
  NOE <- length(events$UTC)                        # Number of Events
  NOI <- NOE - 1                                   # Number of Intervals
  Lambda <- NOE/(events$UTC[NOE] - events$UTC[1])  # Seismic Rate

# and calculate the waiting time between successive events
  Interval <- rep(0, times = NOI)
  for (i in 1:NOI)
    {
    Interval[i] <- (events$UTC[i+1] - events$UTC[i])
    }

# Calculate the density of waiting time intervals
  DensityT <- hist(Interval,breaks=My_Breaks,plot=F)

  Current_Densities <- DensityT$counts/((NOE*Lambda)*(Max_Timespan*2^-seq(15,0)-c(0,Max_Timespan*2^-seq(15,1))))
  Current_Times <- DensityT$mids*Lambda

  Index_L <- 1 + 16*(Counter-1)
  Index_U <- 16*Counter
  Waiting_Data$Time[Index_L:Index_U] <- Current_Times
  Waiting_Data$Density[Index_L:Index_U] <- Current_Densities
  }

Waiting_Data <- subset(Waiting_Data, Density > 0)

#pdf(file="nws_etas_.pdf",paper="a4",width=0,height=0,pointsize=10)
plot(Waiting_Data,log="xy",xlab=expression("T" * lambda),ylab=expression("P" / lambda),main="North West Shelf Earthquakes")

grid(lty=2,col=5)

s <- 10^seq(-4,2,by=0.1)
lines(s,ETAS(s,0.75,0.9,0.04,10),col="red")
lines(s,Corral(s,1.10,1.3,0.75),col="blue")

#dev.off()

# Fit the Corral modified Gamma distribution
#Corral_Fit <- nls(Density ~ Corral(Time, C, alpha, g), data = Waiting_Data, start = list(C=1.25, alpha = 1.5, g = 0.66))


# Fit the ETAS function to the data
#ETAS_Fit <- nls(Density ~ ETAS(Time, alpha, n, theta, rho), data = Waiting_Data, start = list(alpha = 0.66, n = 1.5, theta = 0.04, rho = 100),nls.control(warnOnly=T),trace=TRUE)


