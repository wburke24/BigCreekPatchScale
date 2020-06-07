# modify clim data

setwdhere()
options(stringsAsFactors = FALSE)

library(RHESSysIOinR)

# ========== fix clim data issues ==========
new_clim = read_clim("clim/source/upGGnew")
# summary(new_clim$date)
# "1942-01-01 00:00:00" - "2016-09-02 00:00:00" 
# tdif = as.numeric(difftime(max(new_clim$date), min(new_clim$date), units = "weeks"))/52.25
# 74.6 years

tdif = new_clim$tmax - new_clim$tmin
tbad = which(tdif < 0)
new_clim[tbad,]

tmaxw = new_clim$tmax
tminw = new_clim$tmin

# adjust based on neighboring values
# day 2227 - just barely different
# new_clim[2220:2230,]
tmaxw[tbad[1]] = -3.000000
# day 19604 - looks switched
# new_clim[19600:19608,]
tmaxw[tbad[2]] = new_clim$tmin[tbad[2]]
tminw[tbad[2]] = new_clim$tmax[tbad[2]]
# day 21178 - added negative to tmax
# new_clim[21172:21180,]
tmaxw[tbad[3]] = abs(tmaxw[tbad[3]])

rbad = which(new_clim$rain < 0)
# rain looks fine, as far as data goes
pcpw = new_clim$rain

# min(new_clim$date)
writedate = "1942 01 01 01"

writeLines(text = c(writedate, pcpw), con = "clim/upGG_42_16.rain")
writeLines(text = c(writedate, tmaxw), con = "clim/upGG_42_16.tmax")
writeLines(text = c(writedate, tminw), con = "clim/upGG_42_16.tmin")


# ========== Pad clim (for throwout years for early run) ==========

clim = read_clim("clim/upGG_42_16")

# trim so we only have complete water years - wy 1943 thru 2015
clim = clim[clim$date >= "1942-10-01" & clim$date < "2015-10-01",]

# repeat 1942-10-01 -> 1945-10-01 , add to front
clim_ext = rbind(clim[clim$date >= "1942-10-01" & clim$date < "1945-10-01",], clim )

writedate = "1939 10 01 01" 
writeLines(text = c(writedate, clim_ext$rain), con = "clim/upGG_39_15.rain")
writeLines(text = c(writedate, clim_ext$tmax), con = "clim/upGG_39_15.tmax")
writeLines(text = c(writedate, clim_ext$tmin), con = "clim/upGG_39_15.tmin")
#copy the base file too



# adjust the climate data + 2c
clim_data = "clim/upGG_39_15"

tmin = readLines(paste0(clim_data,".tmin"))
tmax = readLines(paste0(clim_data,".tmax"))
tmin_adj = c(tmin[1], as.character(as.numeric(tmin[2:length(tmin)]) + 2))
tmax_adj = c(tmax[1], as.character(as.numeric(tmax[2:length(tmax)]) + 2))
writeLines(text = tmin_adj,con = paste0(clim_data,"_2C.tmin"))
writeLines(text = tmax_adj,con = paste0(clim_data,"_2C.tmax"))
file.copy(from = paste0(clim_data,".rain"), to = paste0(clim_data,"_2C.rain"), overwrite = TRUE)

file.copy(from = paste0(clim_data,".base"), to = paste0(clim_data,"_2C.base"), overwrite = TRUE)
# replace old name


