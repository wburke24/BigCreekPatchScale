# scenarios
setwdhere()

#devtools::install_github("RHESSys/RHESSysPreprocessing")
#devtools::install("../rhessys/RHESSysPreprocessing")
library(RHESSysPreprocessing)
library(data.table)

# -------------------- Setup Options --------------------


# ---------- Make Scenarios Table ----------

#clim_opts = c("baseline", "clim_change")
clim_opts = c("clim/upGG_39_15", "clim/upGG_39_15_2C")

treatments = c(
  "no treatment",
  "under_0.1 1mo_rm100litter",
  "under_0.25 1mo_rm100litter",
  "under_0.4 1mo_rm100litter",
  "over_0.1 carbon_remains",
  "over_0.25 carbon_remains",
  "over_0.4 carbon_remains",
  "over_0.1 carbon_removed",
  "over_0.25 carbon_removed",
  "over_0.4 carbon_removed",
  "litter_0 prfire"
)
treat_type = c(NA, "under", "under", "under", "over", "over", "over", "over", "over", "over", NA)
#treat_type = c(NA, "shrub", "shrub", "shrub", "conifer", "conifer", "conifer", "conifer", "conifer", "conifer", NA)
treat_intensity = c(NA, 0.1, 0.25, 0.4, 0.1, 0.25, 0.4, 0.1, 0.25, 0.4, NA)
presc_fire = c(F, T, T, T, F, F, F, F, F, F, T)
tec_type = c(NA, "thin_remain", "thin_remain", "thin_remain", "thin_remain", "thin_remain", "thin_remain", "thin_harvest", "thin_harvest", "thin_harvest", NA)
treat = data.frame(treatments, treat_type, treat_intensity, presc_fire, tec_type, stringsAsFactors = FALSE)

#interval = c("1 time", "5 year", "10 year")
interval = c("5", "10", "30")



# min30 max30 med30
# midpoint  2000  1968  1957
# start_wy  1986  1954  1943
# end_wy    2015  1983  1972
# 3 year padding for spin

start_period = c("wet","variable","dry")
# start_date = c("1977 10 01 01", "1943 10 01 01", "2007 10 01 01")
# end_date = c("1997 10 01 01", "1963 10 01 01", "2027 10 01 01")

start_date = c("1953 10 01 01", "1942 10 01 01", "1985 10 01 01")
end_date = c("1983 10 01 01", "1972 9 30 24", "2015 10 01 01")

#end_date = c("1984 9 30 24", "1973 9 30 24", "2016 9 30 24")

# c("wet", as.numeric(difftime(as.Date("1983-09-30"), as.Date("1953-10-01"), unit="days")))
# c("variable", as.numeric(difftime(as.Date("1972-09-30"), as.Date("1942-10-01"), unit="days")))
# c("dry", as.numeric(difftime(as.Date("2015-09-30"), as.Date("1985-10-01"), unit="days")))

# length(seq.Date(as.Date("1953-10-01"), as.Date("1983-09-29"), "day")) *4
# c("variable", as.numeric(difftime(as.Date("1972-09-30"), as.Date("1942-10-01"), unit="days")))
# c("dry", as.numeric(difftime(as.Date("2015-09-30"), as.Date("1985-10-01"), unit="days")))

period = data.frame(start_period, start_date, end_date, stringsAsFactors = FALSE)

soils = c("shallow", "medium", "deep")

sharing = c(0,0.25,0.5,0.75,1)

# worldfiles = c("worldfiles/BC_shrub_N.world", "worldfiles/BC_shrub_S.world", 
#                "worldfiles/BC_shrub_conifer_N.world", "worldfiles/BC_shrub_conifer_S.world", 
#                "worldfiles/BC_conifer_N.world", "worldfiles/BC_conifer_S.world")

veg = c("shrub", "conifer", "shrub_conifer")
aspect = c("N", "S")

scenarios = expand.grid(veg, aspect, clim_opts, treatments, interval, start_period, soils, sharing, stringsAsFactors = FALSE)
names(scenarios) = c("veg", "aspect", "climate", "treatments", "interval", "start_period", "soils", "sharing")

# joins to add data
scenarios = merge(scenarios, treat, by = "treatments", no.dups = FALSE, sort = FALSE)
scenarios = merge(scenarios, period, by = "start_period",  no.dups = FALSE, sort = FALSE)

# ----- clean up invalid runs -----
# get rid of no treatment + 5 and 10 year intervals
scenarios = scenarios[!(scenarios$treatments == "no treatment" & scenarios$interval %in% c("5", "10")  ),]
scenarios[scenarios$treatments == "no treatment",]$interval = "none"

# get rid of overstory thinning + shrub only
scenarios = scenarios[!(scenarios$veg == "shrub" & (!is.na(scenarios$treat_type) & scenarios$treat_type == "over")), ]

# sort by veg
scenarios = scenarios[with(scenarios, order(veg, treatments, soils, climate, start_period, sharing, interval, aspect)), ]

# redo rownums
rownames(scenarios) = 1:length(scenarios[,1])

# add worldfile
scenarios$worldfile = paste0("worldfiles/veg_soils_spin/BC_", scenarios$veg, "_", scenarios$aspect, "_", scenarios$soils, "_", scenarios$sharing, ".world")

save(scenarios, file = "scenarios.rdata")

# ----- shorten scenarios -----
setwdhere()
load("scenarios.rdata")

scen = as.data.table(scenarios)

# start_period: wet, dry, variable
# treatments: method+intensity, presc fire, no treat
# interval chr: 5, 10, 20, "none"
# soils: shallow, medium, deep
# sharing: 0, 0.25, 0.5, 0.75, 1
# treat_type: under, over, none
scen[is.na(scen$treat_type), ]$treat_type = "none"
# treat_intensity: 0.1, 0.25, 0.4, none
scen$treat_intensity =  as.character(scen$treat_intensity)
scen[is.na(scen$treat_intensity), ]$treat_intensity = "none"
# presc_fire: T/F
# veg: conifer, shrub, shrub_conifer
# aspect: N, S
# clim_chg: T/F
scen$clim_chg = scenarios$climate == "clim/upGGmod_2C"
# run: 1:n
scen$run = c(1:nrow(scen))

scen[, c("worldfile", "climate", "tec_type", "start_date", "end_date") ] = NULL

scen$interval =  as.character(scen$interval)

scen = data.table(scen[,"run"],as.data.table(lapply(scen[,!"run"], factor)))

save(scen, file = "scen.rdata")




