# clim

setwdhere()

library(RHESSysIOinR)
library(ggplot2)
library(zoo)
library(plotly)
library(tidyverse)

# ggplot things to try
# use the builder in addins menu
library(gganimate)
library(ggthemes)
library(cowplot)
library(esquisse)

options(stringsAsFactors = FALSE)

clim = read_clim("clim/upGG_42_16")

# plots

dates = data.frame(clim = c(min(clim$date), max(clim$date)),
                   row.names = c("min", "max"))

# trim so we only have complete water years - wy 1943 thru 2015
clim = clim[clim$date >= "1942-10-01" & clim$date < "2015-10-01",]

dates = data.frame(clim = c(min(clim$date), max(clim$date)),
                   row.names = c("min", "max"))


# reminder that there's NAs
any(is.na(clim$tmin) | is.na(clim$tmax) | is.na(clim$rain))


clim$tavg = (clim$tmax + clim$tmin) / 2

clim_yr = clim %>%
  group_by(wy) %>%
  summarise(rain = sum(rain), tmin = mean(tmin), tmax = mean(tmax), tavg = mean(tavg))

mean(clim_yr$rain)*1000
mean(clim$tavg)




clim_yr$rollavg10 = rollmean(clim_yr$rain, 10, fill = NA)
clim_yr$rollavg30 = rollmean(clim_yr$rain, 30, fill = NA)

clim_yr$avg10rank = rank(clim_yr$rollavg10,na.last = "keep")
clim_yr$avg30rank = rank(clim_yr$rollavg30,na.last = "keep")


theme_set(theme_minimal())

p = ggplot(clim_yr, aes(wy)) + 
  #geom_line(aes(y = rain, colour = "rain")) + 
  geom_col(aes(y = rain, fill = avg30rank)) +
  geom_line(aes(y = rollavg10, colour = "rolling 10 year avg"), size = 2) +
  geom_line(aes(y = rollavg30, colour = "rolling 30 year avg"), size = 2) +
  geom_abline(intercept =  median(clim_yr$rain, na.rm = TRUE), slope = 0)

ggplotly(p)


yr_choice = data.frame(min30 = clim_yr$wy[which.min(clim_yr$rollavg30)])
yr_choice$min30 = clim_yr$wy[which.min(clim_yr$rollavg30)]
yr_choice$max30 = clim_yr$wy[which.max(clim_yr$rollavg30)]
yr_choice$med30 = clim_yr$wy[which.min(abs(clim_yr$rollavg30 - median(clim_yr$rollavg30, na.rm = TRUE)))]

#values
min(clim_yr$rollavg30,na.rm = T) * 1000
max(clim_yr$rollavg30, na.rm = T) * 1000
median(clim_yr$rollavg30, na.rm = TRUE) * 1000
mean(clim_yr$rollavg30,na.rm = T)

# range: +15 from yr
yr_choice[2,] = yr_choice[1,] - 14
yr_choice[3,] = yr_choice[1,] + 15

rownames(yr_choice) = c("midpoint", "start_wy", "end_wy")

# med = rank 22  = 1957
# rank 23 = 1975
# rank 21 = 1994

#yr_choice$med30alt = c(1975, NA, NA)



