################################################################################
# 
# Calculate the projected earthquake magnitudes with the Pisarenko-Sornette 
# approach of Extreme Value Theory
#
# Stefan Revets
#
# 2012-01-20 Version 1.0
# 2012-01-23 Version 1.1  Incorporated bootstrapping procedure
# 2012-01-24              Modified & streamlined datafile input
# 2012-01-30              Added & Improved output graphs
# 2012-02-10              Incorporated Kolmogorov-Smirnov to test Poissonian
#                         nature of the processed data set
# 2012-02-13              Included a magnitude perturbation possibility, to test
#                         sensitivity to catalogue magnitude error
#
################################################################################

# First, read in the data from file. Expected format is:
# UTC(YYYY-MM-DD hh:mm:ss[.ss])  Magnitude Latitude  Longitude Depth
# preferably in ascending UTC order
events <- read.table("earthquakes.dat", header = T, sep = "\t")

## Testing sensitivity to magnitude estimate errors
## Just perturb the actual magnitudes by e(0.0, 0.2)
#My_perturb <- rnorm(length(events$Magnitude), mean=0.0, sd = 0.2)
#events$Magnitude <- events$Magnitude + My_perturb

# De-clustering of the catalogue might improve the Poissonian nature of the
# catalogue data, by eliminating any after-shocks.
# Define two functions used to decide whether an event falls into the
# aftershock window
# first function tests for inclusion in the time window
Time_Check <- function(index_main,index_test)
  {
  T_1 <- events$UTC[index_main]
  T_2 <- events$UTC[index_test]
  if (T_2-T_1 < 10^(-0.31+0.46*events$Magnitude[index_main])) Result <- TRUE else Result <- FALSE
  Result
  }
# second function tests for inclusion in the space window  
Distance_Check <- function(index_main,index_test)
  {
   Lat_1 <- events$Latitude[index_main] * pi/180.0
   Long_1 <- events$Longitude[index_main] * pi/180.0
   Lat_2 <- events$Latitude[index_test] * pi/180.0   
   Long_2 <- events$Longitude[index_test] * pi/180.0
   Threshold <- 10^(-0.85+0.46*events$Magnitude[index_main])
   Distance <- 6371 * acos(sin(Lat_1)*sin(Lat_2) + cos(Lat_1)*cos(Lat_2)*cos(Long_2-Long_1))
   if (Distance < Threshold) Result <- TRUE else Result <- FALSE
   Result
  }

# add a logical field, to flag if the entry will be removed as an aftershock
events <- data.frame(events, Keep=rep(T, times = length(events$UTC)))
  
# Reset the time values to number of days since the first event in the database
events$UTC <- as.numeric(difftime(events$UTC,events$UTC[1], units="days"))
events$UTC <- events$UTC - events$UTC[1]

# Proposed set of main-shock magnitudes
Max_Mag <- 4.5
  
# Step through all the events: check if the magnitude falls in the main shock 
# series, then check the following events to see if these fall in the critical
# time-space window

Index <- 1
while (Index < length(events$UTC)-1)     # Stepping through the data matrix
  {
  Step <- 1                            # Starting the aftershock check sequence
  while (Index+Step < length(events$UTC) && Time_Check(Index,Index+Step)==TRUE && (events$Magnitude[Index] > events$Magnitude[Index+Step]))
    {
    if (Distance_Check(Index,Index+Step)==TRUE && events$Magnitude[Index] > Max_Mag) events$Keep[Index+Step] <- FALSE
    Step <- Step + 1
    }
  Index <- Index + 1
  }

# The de-clustered database of events is the subset of Keep==TRUE
events <- subset(events, Keep==TRUE, select = c(UTC,Magnitude,Longitude,Latitude,Depth))

# Reset the time values to number of days since the first event in the database
delta_T_Max <- events$UTC[length(events$UTC)]-events$UTC[1]
events$UTC <- events$UTC - events$UTC[1]

       
# Calculate the Gutenberg-Richter relationship between Magnitude and numbers
# First, set up the magnitude steps, here between 2.0 and 7.0
GR_M_Range <- seq(2.5, 7.0, by=0.1)

# Set up a dataframe for the results                                                                
GR <- data.frame(Threshold=GR_M_Range,Number=rep(0, times=length(GR_M_Range)))

# Calculate the number of events exceeding the defined magnitude range
for (i in 1:length(GR_M_Range)) 
  {
  GR$Number[i] <- length(subset(events$Magnitude, events$Magnitude > GR_M_Range[i]))
  }

################################################################################
#
# EVD calculations
#

# Ensure the necessary libraries are available
require(evd)

# First, set the lower magnitude cut-off
M_Min <- 3.5

# Now define time steps for the Generalised Extreme Value distribution,
# in days
delta_T <- 10
Time_Steps <- seq(20, 300, by=delta_T)

# Set up the number of data shuffles to bootstrap the GEV parameters
# and improve their accuracy (i.e., reduce variability)
Bootstrap_Total <- 100
shuffle_events <- events
 
# The fitted parameters go into an 3-d array
GEV_Parameters <- array(0, c(length(Time_Steps),4,Bootstrap_Total))

for (Re_runs in 1:Bootstrap_Total)
  {# The idea behind the bootstap approach is that shuffling the magnitudes
  # around amounts to a resampling of the population of the events, whilst
  # maintaining the distribution in time
  shuffle_events$Magnitude <- sample(events$Magnitude)
  # We need to step through the entire events dataset in contiguous blocks
  # of size Time_Steps[i], and determine the maximum magnitude in each of
  # the intervals.

  for (i in 1:length(Time_Steps))
    {# Looping over the Time Intervals
    # Determine the contents of the successive time bins through the hist function:
    # use the number in each bin (and accumulate) to find the position of the data
    # entries in the events matrix
    my_breaks <- seq(0,delta_T_Max+Time_Steps[i],by=Time_Steps[i])
    Time_Hist <- hist(shuffle_events$UTC,breaks=my_breaks,plot=F)
  
    # The numbers in each bin are stored in Time_Hist$counts and can be used now
    # to calculate the maximum magnitude encountered in each of the time bins
    Bin_low <- 0
    Bin_high <- 0
    Bin_Max_Magnitudes <- rep(0, times=length(Time_Hist$counts))
  
    for (Bins in 1:length(Time_Hist$counts))
      {
      Bin_high <- Bin_low + Time_Hist$counts[Bins]
      Bin_Max_Magnitudes[Bins] <- max(shuffle_events$Magnitude[Bin_low:Bin_high])
      if (Bin_Max_Magnitudes[Bins] < M_Min) Bin_Max_Magnitudes[Bins] <- NA
      Bin_low <- Bin_high
      }

    # Calculate the MLE of the GEV distribution and store the results
    # The order is: T, loc, scale, shape, error_loc, error_scale, error_shape
    GEV_Fit <- fgev(Bin_Max_Magnitudes,std.err=F)
    GEV_Parameters[i, ,Re_runs] <- c(Time_Steps[i],fitted.values(GEV_Fit))
    }# End of Time Interval Looping  
  }

# Present the results by performing a statistical summary of the parameter
# estimates

GEV_Results <- array(0, c(length(Time_Steps),10))

for (i in 1:length(Time_Steps))
  { 
  GEV_Results[i,1] <- Time_Steps[i]
  GEV_Results[i,2:4] <- quantile(GEV_Parameters[i,2,],probs=c(0.16,0.50,0.84))
  GEV_Results[i,5:7] <- quantile(GEV_Parameters[i,3,],probs=c(0.16,0.50,0.84))
  GEV_Results[i,8:10] <- quantile(GEV_Parameters[i,4,],probs=c(0.16,0.50,0.84))
  }

Shape <- GEV_Results[,9]  
Scale <- GEV_Results[,6]
Location <- GEV_Results[,3]

# Now we can calculate estimates of maximum magnitudes for arbitrary time
# in the future

Tau <- c(365000*1:5)
Q <- 0.975

Quantiles <- array(0, c(length(Time_Steps),length(Tau)))
for (Tau_i in 1:length(Tau))
  {
  Quantiles[,Tau_i] <- Location + ((Tau[Tau_i]/(log(1/Q)*Time_Steps))^Shape - 1) * Scale / Shape
  }

# Make some plots

# Let's test how well the reduced data set fulfills the Poisson Process
Test_Data <- subset(events, events$Magnitude > M_Min, select=UTC)
Test_Data$UTC  <- Test_Data$UTC - Test_Data$UTC[1]
NOE <- length(Test_Data$UTC)
TL <- Test_Data$UTC[NOE]
# Generate the equivalent Poissonian data set
Poisson_Data <- seq(1, TL, by=TL/NOE)
#  and carry out the Kolmogorov-Smirnov test
Test_Result <- ks.test(events$UTC,Poisson_Data)
       
pdf(file="earthquakes_gev.pdf",paper="a4",width=0,height=0,pointsize=10)
op <- par(mfrow = c(3,2))
plot(GR,log="y",main="Gutenberg-Richter Plot");grid(lty=2,col=5)
abline(v=M_Min,col=2)
My_List <- subset(events, Magnitude > M_Min, select = c(UTC,Magnitude))
plot(My_List$UTC/My_List$UTC[length(My_List$UTC)],type="l",xlab="Event Number",ylab="Normalised Occurrence Time",main="Poisson Fit");grid(lty=2,col=5)
abline(0,1/length(My_List$UTC),col=4)
text(NOE/5,0.8,"p-value:")
text(NOE/5,0.7,round(Test_Result$p.value,5))
matplot(GEV_Results[,1],GEV_Results[,8:10],type="l",lty=1,col=c(1,2,1),main="GEV Parameter Estimation",xlab="T Window (days)", ylab="Shape Parameter");grid(lty=2,col=5)
matplot(GEV_Results[,1],GEV_Results[,5:7],type="l",lty=1,col=c(1,2,1),main="GEV Parameter Estimation",xlab="T Window (days)", ylab="Scale Parameter");grid(lty=2,col=5)
matplot(GEV_Results[,1],GEV_Results[,2:4],type="l",lty=1,col=c(1,2,1),main="GEV Parameter Estimation",xlab="T Window (days)", ylab="Location Parameter");grid(lty=2,col=5)
matplot(Time_Steps,Quantiles,ylim=c(6,9),type="l",lty=1,col=1,main="GEV Maximum Magnitude Estimation",xlab="T Window (days)", ylab="0.975 Magnitude Quantile");grid(lty=2,col=5)
par(op)
#dev.off()
