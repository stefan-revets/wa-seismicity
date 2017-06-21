
library("data.table")
library("dplyr")
library("evd")
library("geosphere")
library("tidyr")

events <- read.table("../data/earthquakes.csv", header = TRUE,
                     sep = ",", stringsAsFactors = FALSE)
events <- select(events, -Sydney.Date, -Sydney.Time,
                 -Approximate.location, -ORIGIN.ID)
events <- unite(events, UTC, c(UTC.Date, UTC.Time), sep = " ")
events$UTC <- as.POSIXct(events$UTC, tz = "GMT")
events <- arrange(events, UTC)

events <- rename(events,
                 Mag = Magnitude,
                 Long = Longitude,
                 Lat = Latitude,
                 Depth = Depth..km.)
events <- mutate(events, Aftershock = FALSE)

pdf(file="../output/time_mag.pdf")
plot(events$UTC,jitter(events$Mag),
     main = "Earthquakes recorded over time",
     xlab = "Date",
     ylab = "Magnitude",
     pch = 20,
     col = "blue")
grid()

dev.off()

pdf(file="../output/mag_issues.pdf")
GR_pre <- events %>%
    filter(year(UTC) < 1992) %>%
    count(Mag) %>%
    mutate(n = cumsum(n)) %>%
    mutate(n = max(n) - n + 1)

GR_post <- events %>%
    filter(year(UTC) > 1991) %>%
    count(Mag) %>%
    mutate(n = cumsum(n)) %>%
    mutate(n = max(n) - n + 1)

GR_corr<- events %>%
    filter(year(UTC) < 1992) %>%
    mutate(Mag = Mag - 0.5) %>%
    count(Mag) %>%
    mutate(n = cumsum(n)) %>%
    mutate(n = max(n) - n + 1)

plot(GR_pre, log = "y", pch = 19,
     main = "Gutenberg-Richter relation")
points(GR_post)
points(GR_corr, col = "green")
grid()

legend(2, 5,
       legend = c("pre 1992","post 1992","pre 1992 corrected"),
       pch = c(19,21,21),
       col = c("black","black","green"))
dev.off()

events$Mag[(year(events$UTC) < 1992)] <-
    events$Mag[(year(events$UTC) < 1992)] - 0.5

pdf(file="../output/ecdf.pdf")
events <- mutate(events,
                 UTC_d = as.numeric(difftime(UTC,UTC[1],
                                             unit = "days")))
plot(ecdf(events$UTC_d/365.4 + year(events$UTC[1])),
     main = "Cumulative distribution of earthquakes",
     xlab = "Year",
     pch = 20,
     col = "blue")
grid()

dev.off()

pdf(file="../output/ecdf_restricted.pdf")
events <- events %>%
    filter(year(UTC) > 1979, year(UTC) < 2003) %>%
    mutate(UTC_d = as.numeric(difftime(UTC,UTC[1], unit = "days")))

plot(ecdf(events$UTC_d/365.4 + year(events$UTC[1])),
     main = "Cumulative distribution of earthquakes",
     xlab = "Year",
     pch = 20,
     col = "blue")
grid()

dev.off()

NOE <- length(events$UTC_d)
TL <- max(events$UTC_d)
Poisson_Data <- seq(1, TL, by=TL/NOE)
KS_Result <- ks.test(events$UTC_d,Poisson_Data)
KS_Result

Max_Mag <- 5.0
main_shocks <- filter(events, Mag > Max_Mag)
main_events <- mutate(main_shocks,
                      delta_T = 10^(-0.31 + 0.46 * Mag),
                      delta_D = 10^(-0.85 + 0.46 * Mag))

events <- filter(events, Mag <= Max_Mag)

Index <- 1
  while (Index < length(events$UTC)-1)
  {
      main_result <- main_events %>%
          mutate(test_dT = as.numeric(difftime(events$UTC[Index], UTC,
                                               unit = "days")),
                 test_dD = distHaversine(cbind(Long, Lat),
                                         c(events$Long[Index], events$Lat[Index]))/1000
                 ) %>%
          filter(test_dT > 0, test_dT < delta_T, test_dD < delta_D)
      if (nrow(main_result) > 0) {events$Aftershock[Index] = TRUE}
      
    Index <- Index + 1
    }

events <- rbind(events, main_shocks)
events <- arrange(events, UTC)
events <- filter(events, Aftershock == FALSE)

M_Min <- 2.5
delta_T <- 10
Time_Steps <- seq(20, 300, by = delta_T)
delta_T_Max <- last(events$UTC_d) - first(events$UTC_d)

Bootstrap_Total <- 100
shuffle_events <- events

GEV_Parameters <- array(0, c(length(Time_Steps),4,Bootstrap_Total))

for (Re_runs in 1:Bootstrap_Total){
        shuffle_events$Mag <- sample(events$Mag)
        for (i in 1:length(Time_Steps)){
            shuffle_events <- shuffle_events %>%
                mutate(block = UTC_d %/% Time_Steps[i])
            Max_Mags <- shuffle_events %>%
                group_by(block) %>%
                summarize(value = max(Mag))
            GEV_Fit <- fgev(Max_Mags$value,std.err=F)
            GEV_Parameters[i, ,Re_runs] <-
                c(Time_Steps[i],fitted.values(GEV_Fit))
        }
}

GEV_Results <- array(0, c(length(Time_Steps),10))

for (i in 1:length(Time_Steps))
  { 
  GEV_Results[i,1] <- Time_Steps[i]
  GEV_Results[i,2:4] <- quantile(GEV_Parameters[i,2,],
                                 probs=c(0.16,0.50,0.84))
  GEV_Results[i,5:7] <- quantile(GEV_Parameters[i,3,],
                                 probs=c(0.16,0.50,0.84))
  GEV_Results[i,8:10] <- quantile(GEV_Parameters[i,4,],
                                  probs=c(0.16,0.50,0.84))
  }

Shape <- GEV_Results[,9]  
Scale <- GEV_Results[,6]
Location <- GEV_Results[,3]

Tau <- c(365000*1:5)
Q <- 0.975

Quantiles <- array(0, c(length(Time_Steps),length(Tau)))
for (Tau_i in 1:length(Tau))
  {
      Quantiles[,Tau_i] <- Location +
          ((Tau[Tau_i]/(log(1/Q)*Time_Steps))^Shape - 1) * Scale / Shape
  }

GR <- events %>%
    count(Mag) %>%
    mutate(n = cumsum(n)) %>%
    mutate(n = max(n) - n + 1)

pdf(file="../output/earthquakes_gev.pdf")
op <- par(mfrow = c(3,2))

plot(GR,log="y",main="Gutenberg-Richter Plot")
grid(lty=2,col=5)
abline(v=M_Min,col=2)
My_List <- subset(events, Mag > M_Min, select = c(UTC_d,Mag))

plot(My_List$UTC_d/My_List$UTC_d[length(My_List$UTC_d)],
     type="l",
     xlab="Event Number",
     ylab="Normalised Occurrence Time",
     main="Poisson Fit")
grid(lty=2,col=5)
abline(0,1/length(My_List$UTC_d),col=4)
text(NOE/5,0.8,"p-value:")
text(NOE/5,0.7,round(KS_Result$p.value,5))

matplot(GEV_Results[,1],GEV_Results[,8:10],
        type="l",
        lty=1,
        col=c(1,2,1),
        main="GEV Parameter Estimation",
        xlab="T Window (days)",
        ylab="Shape Parameter")
grid(lty=2,col=5)

matplot(GEV_Results[,1],GEV_Results[,5:7],
        type="l",
        lty=1,
        col=c(1,2,1),
        main="GEV Parameter Estimation",
        xlab="T Window (days)",
        ylab="Scale Parameter")
grid(lty=2,col=5)

matplot(GEV_Results[,1],GEV_Results[,2:4],
        type="l",
        lty=1,
        col=c(1,2,1),
        main="GEV Parameter Estimation",
        xlab="T Window (days)",
        ylab="Location Parameter")
grid(lty=2,col=5)

matplot(Time_Steps,Quantiles,ylim=c(5,8),
        type="l",
        lty=1,
        col=1,
        main="GEV Maximum Magnitude Estimation",
        xlab="T Window (days)",
        ylab="0.975 Magnitude Quantile")
grid(lty=2,col=5)

par(op)
dev.off()
