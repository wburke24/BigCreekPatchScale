# collect output

library(data.table)
library(RHESSysIOinR)
#library(notifier)
#library(googledrive)

setwdhere()

load("scenarios.rdata")

veg = unique(scenarios$veg)
folders = paste0("output/", veg)
# get from run script 
vars = output_variables$variable
scenarios$run = 1:nrow(scenarios)


fixdates =F
if (fixdates) {
  simnum = NULL
  dates = NULL
  simday = NULL
  ct = 0
  for (r in scenarios[scenarios$veg == "conifer", ]$run ) {
    ct = ct+1
    dseq = seq.POSIXt(from = as.POSIXct(gsub(" ","-",scenarios$start_date[r])), to = as.POSIXct(gsub(" ","-",scenarios$end_date[r])) - (60*60*24), by = "DSTday")
    dates = c(dates, dseq)
    simnum = c(simnum, rep(r, times = length(dseq)))
    simday = c(simday, c(1:length(dseq)))
    cat(ct, "\n")
  }
  dt_conifer = data.table(simnum, dates, simday)
  
  simnum = NULL
  dates = NULL
  simday = NULL
  ct = 0
  for (r in scenarios[scenarios$veg == "shrub", ]$run ) {
    ct = ct+1
    dseq = seq.POSIXt(from = as.POSIXct(gsub(" ","-",scenarios$start_date[r])), to = as.POSIXct(gsub(" ","-",scenarios$end_date[r])) - (60*60*24), by = "DSTday")
    dates = c(dates, dseq)
    simnum = c(simnum, rep(r, times = length(dseq)))
    simday = c(simday, c(1:length(dseq)))
    cat(ct, "\n")
  }
  dt_shrub = data.table(simnum, dates, simday)
  
  simnum = NULL
  dates = NULL
  simday = NULL
  ct = 0
  for (r in scenarios[scenarios$veg == "shrub_conifer", ]$run ) {
    ct = ct+1
    dseq = seq.POSIXt(from = as.POSIXct(gsub(" ","-",scenarios$start_date[r])), to = as.POSIXct(gsub(" ","-",scenarios$end_date[r])) - (60*60*24), by = "DSTday")
    dates = c(dates, dseq)
    simnum = c(simnum, rep(r, times = length(dseq)))
    simday = c(simday, c(1:length(dseq)))
    cat(ct, "\n")
  }
  dt_shrub_confier = data.table(simnum, dates, simday)
  
  
  # make into datatable
  dt_conifer_rep = as.data.table(lapply(dt_conifer, rep, each = 4))
  dt_shrub_rep = as.data.table(lapply(dt_shrub, rep, each = 4))
  dt_shrub_confier_rep = data.table()
  dt_shrub_confier_rep[ , simrun := rep(dt_shrub_confier$simnum, each = 9)]
  dt_shrub_confier_rep[ , dates := rep(dt_shrub_confier$dates, each = 9)]
  dt_shrub_confier_rep[ , simday := rep(dt_shrub_confier$simday, each = 9)]
  # add patch ID
  dt_conifer_rep[ , patchID := rep(c(1:4), length.out = nrow(dt_conifer_rep))]
  dt_shrub_rep[ , patchID := rep(c(1:4), length.out = nrow(dt_shrub_rep))]
  dt_shrub_confier_rep[ , patchID := rep(c(1:9), length.out = nrow(dt_shrub_confier_rep))]
  
  save(dt_conifer_rep, file = "conifer_date.rdata")
  save(dt_shrub_rep, file = "shrub_date.rdata")
  save(dt_shrub_confier_rep, file = "shrub_confier_date.rdata")
  #rename to dt_date
}
# 
# load("conifer_date.rdata")
# dt_date = dt_conifer_rep
# save(dt_date, file = "conifer_date.rdata")
# 
# load("shrub_date.rdata")
# dt_date = dt_shrub_rep
# save(dt_date, file = "shrub_date.rdata")
# 
# load("shrub_confier_date.rdata")
# dt_date = dt_shrub_confier_rep
# save(dt_date, file = "shrub_confier_date.rdata")

# # undo the wide
# pnum = 4
# if (folders[j] == "output/shrub_conifer") {pnum = 9}
# 
# origlen = (pnum * 10956 * (ncol(read_in)/3) ) + (pnum * 10957 * (ncol(read_in)/3)) + (pnum * 10956 * (ncol(read_in)/3))
# # badlen =  dim(read_in)[1] * dim(read_in)[2]
# read_longbad = matrix(as.matrix(read_in), ncol = 1)
# read_long = read_longbad[1:origlen]
# 
# rm(read_in, read_longbad)
# gc()
# 
# load(paste0(substr(folders[j],8,nchar(folders[j])), "_date.rdata" ))
# 
# if(nrow(dt_date) != length(read_long)) {
#   dt_date[, lapply()]
#   #a[, paste0("s", cols.to.sum) := lapply(.SD, sum), by = id, .SDcols = cols.to.sum]
# }
# 
# # combine here
# dt_date[ , value := read_long]
# 
# rm(read_long)
# gc()
# 
# runs_long = scenarios$run[scenarios$start_period == "variable" & scenarios$run %in% unique(dt_date$simnum)]
# dt_date[ simnum %in% runs_long & simday == 366, ]

vars = vars[6:11]

for (i in vars) {
  
  DT_all = NULL
  
  for (j in seq_along(folders)) {
    
    # get the data
    read_in = fread(file.path(folders[j], i)) 
    
    # set colnames to just the nums
    colnames(read_in)[startsWith(colnames(read_in), "run")] = which(scenarios$veg == veg[j])
    
    # add dates to data
    dt = add_dates(read_in)
    rm(read_in)
    gc()
    
    # get full output with areas
    allout_read = readin_rhessys_output(folders[j])
    
    # merge to get areas
    if (!"area" %in% colnames(dt)) {
      dt = merge.data.table(dt, allout_read$pd[,c("date", "patchID", "area")], by = c("date", "patchID"))
    }

    dt_agg = patch_fam_agg(X = dt, na.rm = TRUE) 
    rm(dt)
    gc()
        
    dt_melt = melt(data = dt_agg, measure.vars = which(!colnames(dt_agg) %in% c("familyID", "date", "stratumID")), variable.name = "run", variable.factor = FALSE)
    dt_melt[, run := as.integer(run), ]
    #dt_melt$run = as.integer(dt_melt$run)
    #colnames(dt_melt)[colnames(dt_melt) == "value"] = i
    rm(dt_agg)
    gc()
    
    if (j == 1) {
      DT_all = dt_melt
    } else {
      DT_all = rbind(DT_all, dt_melt)
    }
    rm(dt_melt)
    gc()
    
  }
  
  save(DT_all, file = paste0("output/R_obj/", i, ".rdata"))
  
}

beepr::beep(sound = 2)

notifier::notify(
  title = "R",
  msg = c("Code is done running")
)

#system('7z e -o <output_dir> <archive_name>')

# subset outputs for HydroShare to get <1GB per object
sub_out = F
if (sub_out) {
  
  # pct_chg, scenarios, are already small enough
  # daily npp, standc, et, pspread, and heights (seperate for shrub and conifer+shrub)
  load("output/R_obj/psn.rdata")
  npp = DT_all[,c(3,4)]
  save(npp, file = "output/HydroShare_out/npp.rdata", compress = "xz")
  rm(DT_all, npp)
  gc()
  
  load("output/R_obj/fire_et.rdata")
  et = DT_all[,c(3,4)]
  save(et, file = "output/HydroShare_out/et.rdata", compress = "xz")
  rm(DT_all, et)
  gc()
  
  load("output/R_obj/plantc.rdata")
  standc = DT_all[,c(3,4)]
  save(standc, file = "output/HydroShare_out/standc.rdata", compress = "xz")
  rm(DT_all)
  gc()
  
  load("output/R_obj/pspread.rdata")
  firespread = DT_all[,c("run","value")]
  save(firespread, file = "output/HydroShare_out/firespread.rdata", compress = "xz")
  rm(DT_all)
  gc()
  
  load("output/R_obj/height.rdata")
  load("scen.rdata")
  
  heights = merge.data.table(DT_all, scen[,c("veg", "run")], by = "run", all = FALSE)
  shrub_heights = heights[heights$veg == "shrub", c("run", "value")]
  conifer_shrub_heights = heights[heights$veg == "conifer", c("run", "date", "stratumID", "value")]

  rm(DT_all, heights)
  gc()
  
  save(shrub_heights, file = "output/HydroShare_out/shrub_heights.rdata", compress = "xz")
  save(conifer_shrub_heights, file = "output/HydroShare_out/conifer_shrub_heights.rdata", compress = "xz")
  
  rm(conifer_shrub_heights, shrub_heights)
  gc()
  
  # test
  load("output/HydroShare_out/npp.rdata")
  load("output/HydroShare_out/et.rdata")
  load("output/HydroShare_out/standc.rdata")
  
  
  
  
}



