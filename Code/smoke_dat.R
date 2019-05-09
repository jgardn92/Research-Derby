install.packages("lubridate")
install.packages("chron")

require(lubridate)
require(chron)

smoke_days <- read.csv("Data/smoke_template.csv",header = TRUE)
dates <- as.Date(as.POSIXct(smoke_days$Date,format = "%m/%d/%Y"))

julian_day <- yday(dates)
year <- year(dates)

smoke_days2 <- cbind(smoke_days,julian_day, year)
smoke_days2 <- smoke_days2[,2:5]

saveRDS(smoke_days2,"Output/smoke_days.RDS")
