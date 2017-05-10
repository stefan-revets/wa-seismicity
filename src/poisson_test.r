################################################################################
#
# Test if a particular time series is Poissonian or not, using
# Kolmogorov-Smirnov
#
# 2012-02-10 Stefan Revets, version 1.0
#
################################################################################

# Expect the following variables
#
# Test_Data <- (Time, Value)
#   or Test_Data <- data.frame(Time=events$UTC,Value=events$Magnitude)
# NOE <- number of events  (=length(Test_Data$Time))
# TL <- length of entire time interval (=Test_Data$Time[NOE} - Test_Data$Time[1])

# Generate a Poissonian data set
Poisson_Data <- seq(1, TL, by=TL/NOE)

# Carry out the Kolmogorov-Smirnov test
Test_Result <- ks.test(Test_Data$Time,Poisson_Data)
qqplot(Test_Data$Time,Poisson_Data,main=Test_Result$method,xlab="Test Data",ylab="Poissonian Data",type="l",col=2)
abline(0,1,col=4)
