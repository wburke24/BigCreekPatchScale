setwdhere()

library(RHESSysIOinR)
library(RHESSysPreprocessing)
library(data.table)
library(tictoc)

options(stringsAsFactors = FALSE) # just in case


# -------------------- //// New spin up from initialized soil stores (spin_soil_CN) //// --------------------

# spin separately for each veg, aspect, soil depth, sharing coef
scenarios = NULL
load("scenarios.rdata")

veg = unique(scenarios$veg)
aspect = unique(scenarios$aspect)
soils  = unique(scenarios$soils)
sharing = unique(scenarios$sharing)

spin_scenarios = expand.grid(veg, aspect, soils, sharing, stringsAsFactors = FALSE)
names(spin_scenarios) = c("veg", "aspect", "soils", "sharing")
spin_scenarios$worldfile_out = paste0("worldfiles/veg_soils_spin/BC_", spin_scenarios$veg, "_", spin_scenarios$aspect, "_", spin_scenarios$soils, "_", spin_scenarios$sharing, ".world")
spin_scenarios$worldfile_in = paste0("worldfiles/spin_soil_CN/BC_", spin_scenarios$veg, "_", spin_scenarios$aspect, ".world")

spin_scenarios$run = c(1:nrow(spin_scenarios))
spin_scenarios = as.data.table(spin_scenarios)

save(spin_scenarios, file = "spin_scenarios.rdata")

# -------------------- Start runs --------------------

# -------------------- Computer ID --------------------
comp_ID = Sys.info()["nodename"]

# -------------------- Multiscale dev stuff --------------------
compile = FALSE
if (compile) {
  comp_out = compile_rhessys(location = "../../../rhessys/rhessys-develop/", delete_objs = TRUE, 
                             destination = file.path("bin",comp_ID), CFLAGS = "-w -O2")
}

# -------------------- Project/Run Name --------------------
name = "BC_patch"
n_sim = 1

# -------------------- Parameter Method --------------------
parameter_method = "exact_values"

# spin try 1
#start_end = c("1942 10 01 01", "2002 10 01 24")

# spin try 2
clim_data = "clim/upGGmod"
clim_in = clim_data
clim_out = "clim/upGGmod_repeated"
n_reps = 5

clim_repeat = function(clim_in, clim_out, n_reps) {
  # using this so clim could be repeated using dates in the future
  read = read_clim(clim_in = clim_in)
  clim_rep = do.call("rbind", replicate(n_reps, read, simplify = FALSE))
  cols_out = colnames(clim_rep)[!colnames(clim_rep) %in% c("date", "year", "month", "day", "wy", "yd", "wyd") ]
  files_out = paste0(clim_out, ".", cols_out)
  
  # need an IF in case lengths are mismatches
  date_head = paste(clim_rep$year[1], clim_rep$month[1], clim_rep$day[1], "1")
  write_cols = rbind(date_head, clim_rep[,cols_out] )
  shh = mapply(write, write_cols, files_out)
  
  # copy edit rename basestation
  read_base = readLines(paste0(clim_in, ".base"))
  read_base[endsWith(read_base, "\tclimate_prefix")] = gsub(clim_in, clim_out, read_base[endsWith(read_base, "\tclimate_prefix")])
  writeLines(read_base, paste0(clim_out, ".base"))
}  

clim_repeat(clim_in, clim_out, n_reps)
clim_data = clim_out
start_end = c("1942 10 01 01", "2092 10 01 24")

# -------------------- Input RHESSys --------------------
input_rhessys <- list()
input_rhessys$rhessys_version <- file.path("bin",comp_ID,"rhessys7.1.1")
input_rhessys$tec_file <- "tecfiles/BC.tec"
input_rhessys$world_file <- NULL
input_rhessys$world_hdr_prefix <- "BC"
input_rhessys$flow_file <- NULL
input_rhessys$start_date <- start_end[1]
input_rhessys$end_date <- start_end[2]
input_rhessys$output_folder <- "output/spinup/"
input_rhessys$output_filename <- "BC"
input_rhessys$command_options <- c("-g -p -c -msr") #-v -6  

# -------------------- Input Headers --------------------
input_hdr_list <- list()
input_hdr_list$basin_def <- c("defs/basin_p301.def")
input_hdr_list$hillslope_def <- c("defs/hill_p301.def")
input_hdr_list$zone_def <- c("defs/zone_p301.def")
input_hdr_list$soil_def <- c("defs/soil_forestshrub.def")
input_hdr_list$landuse_def <- c("defs/lu_p301.def")
input_hdr_list$patch_def <- c("defs/patch_p301.def")
input_hdr_list$stratum_def <- c("defs/veg_p301_shrub.def", "defs/veg_p301_conifer.def")
input_hdr_list$base_stations <- c("clim/Grove_lowprov_clim.base")

# -------------------- Parameters --------------------
input_preexisting_table <- NULL

# -------------------- Soils --------------------
soilslist = NULL
load("soilslist.rdata")

# -------------------- Standard (soil/subsurface) Parameters--------------------
input_standard_par_list <- list(m = c(2), k = c(2), m_v = c(2), k_v = c(2), pa = c(1.15), po = c(0.766), gw1 = c(0.24), gw2 = c(0.2))

# -------------------- Make Climate Basestation --------------------
#input_clim_base_list <- NULL

input_clim_base_list = clim_auto(base_station_id = 101, x_coordinate = 100.0, y_coordinate = 100.0, z_coordinate = 1748, 
                                 effective_lai = 3.5, screen_height = 2, daily_prefix = clim_data)

# -------------------- Make Dated Sequence --------------------
input_dated_seq_list <- NULL

# -------------------- Make Tec File --------------------
# tec handled inside iteration loop
#start_end = read_clim(clim_data, dates_out = TRUE)
#start_end = c("1942 10 01 01", "2002 10 01 01")
input_tec_data = input_tec(NULL, start_end = start_end)

# -------------------- Output --------------------
output_method = "r"

# patch/canopy level vars
output_variables <- data.frame(out_file = character(), variable = character(), stringsAsFactors = FALSE)
output_variables[1,] <- data.frame("pd", "Qout")
output_variables[2,] <- data.frame("pd", "psn")
output_variables[3,] <- data.frame("pdg", "lai")
output_variables[4,] <- data.frame("pdg", "plantc")
output_variables[5,] <- data.frame("cd", "height")
output_variables[6,] <- data.frame("cdg", "leafc")
output_variables[7,] <- data.frame("pd", "fire_et")
output_variables[8,] <- data.frame("pd", "pet")

# ========================= Iterate Scenarios =========================


# ||||| -- RUN FROM HERE -- |||||

tic()

veg = unique(spin_scenarios$veg)

for (v in seq_along(veg)) {
  
  runs = which(spin_scenarios$veg == veg[v])
  
  input_rhessys$output_filename = veg[v]
  
  if (!dir.exists("output/spinup/allsim")) {
    dir.create("output/spinup/allsim")
  }
  
  # run count
  runct = 0
  
  for (i in runs) {
    
    runct = runct + 1
    cat("\n ========== Scenario ",runct," ========== \n\n")
    
    # ----- input_rhessys vars -----
    input_rhessys$world_file = spin_scenarios$worldfile_in[i]
    input_rhessys$flow_file = paste0("flowtables/BC_", spin_scenarios$veg[i], ".flow")
    
    # ----- veg, soil and sharing parameters -----
    input_def_list = c(list(list(input_hdr_list$landuse_def[1], "sh_l", spin_scenarios$sharing[i]),
                            list(input_hdr_list$stratum_def[1], "epc.height_to_stem_coef", 2.0),
                            list(input_hdr_list$stratum_def[2], "epc.height_to_stem_coef", 6.0)),
                       soilslist[[ spin_scenarios$soils[i] ]] )
    
    # -------------------- End Run Setup--------------------
    # -------------------- Run RHESSys --------------------
    run_rhessys(parameter_method = parameter_method,
                output_method = output_method,
                input_rhessys = input_rhessys,
                input_hdr_list = input_hdr_list,
                input_preexisting_table = input_preexisting_table,
                input_def_list = input_def_list,
                input_standard_par_list = input_standard_par_list,
                input_clim_base_list = input_clim_base_list,
                input_dated_seq_list = input_dated_seq_list,
                input_tec_data = input_tec_data,
                output_variables = NULL, 
                return_data = FALSE)
    
    # -------------------- Get Output --------------------
    
    #substr(scenarios$worldfile[i], 12, nchar(scenarios$worldfile[i]) - 6)
    
    data_out = select_output_variables_R(output_variables = output_variables,
                                         output_folder = input_rhessys$output_folder,
                                         output_filename = input_rhessys$output_filename,
                                         run = runct,
                                         max_run =  length(runs),
                                         return_data = FALSE)
    
    # copy the output current state worldfile
    # file.rename(from = paste0("worldfiles/spin_soil_CN/BC_", spin_scenarios$veg[i], "_", spin_scenarios$aspect[i], ".world.Y2002M10D1H1.state"), 
    #           to = spin_scenarios$worldfile_out[i])
    
    file.rename(from = paste0("worldfiles/spin_soil_CN/BC_", spin_scenarios$veg[i], "_", spin_scenarios$aspect[i], ".world.Y2092M10D1H1.state"), 
                to = spin_scenarios$worldfile_out[i])
    
    cat("\n ========== End Scenario ",runct," ========== \n\n")
    
  }
  
  file.rename("output/spinup/allsim", paste0("output/spinup/", veg[v]))
  
  
}

toc()

# ==================== end spinup scenario iterations ====================
# 

# spinup output collection


folders = paste0("output/spinup/", veg)
#rh_out_pre = c("output/spinup/BC_conifer", "output/spinup/BC_shrub_conifer", "output/spinup/BC_shrub")

vars = output_variables$variable

for (i in vars) {
  
  DT_all = NULL
  
  for (j in seq_along(folders)) {
    
    # get the data
    read_in = fread(file.path(folders[j], i)) 
    
    # set colnames to just the nums
    colnames(read_in)[startsWith(colnames(read_in), "run")] = which(spin_scenarios$veg == veg[j])
    
    # add dates to data
    dt = add_dates(read_in)
    
    # get full output with areas
    allout_read = readin_rhessys_output(folders[j])

    # merge to get areas
    if (!"area" %in% colnames(dt)) {
      dt = merge.data.table(dt, allout_read$pd[,c("date", "patchID", "area")], by = c("date", "patchID"))
    }
    
    dt_agg = patch_fam_agg(X = dt, na.rm = TRUE) 
    
    dt_melt = melt(data = dt_agg, measure.vars = which(!colnames(dt_agg) %in% c("familyID", "date", "stratumID")), variable.name = "run", variable.factor = FALSE)
    dt_melt$run = as.integer(dt_melt$run)
    #colnames(dt_melt)[colnames(dt_melt) == "value"] = i
    
    if (j == 1) {
      DT_all = dt_melt
    } else {
      DT_all = rbind(DT_all, dt_melt)
    }
    
  }
  
  save(DT_all, file = paste0("output/spinup/R_obj/", i, ".rdata"))
  
  rm(read_in, dt, dt_agg, dt_melt)
  gc()
  
}



# plot the outputs

# ----- load data, basic processing -----

library(tidyverse)
library(gridExtra)
# read data and add dates for aggregation

read_robj = function(path) {
  DT_all = NULL
  load(path)
  if ("familyID" %in% colnames(DT_all)) {
    DT_all$familyID = NULL
  }
  days = length(DT_all$run[DT_all$run == 1])
  DT_all[, sim_day :=  rep(c(1:days), nrow(DT_all)/days), ]
  DT_all[, mn := lubridate::month(date), ]
  DT_all[, yr := lubridate::year(date), ]
  DT_all[, wy := data.table::fifelse(mn >= 10, yr + 1, yr), ]
  DT_all[, wym := data.table::fifelse(mn >= 10, mn - 9, mn + 3), ]
  DT_all[, sim_yr := wy - min(DT_all$yr), ]
  #beepr::beep(2)
  return(DT_all)
}


# MONTHLY
mean2mn = function(DT, scenarios) {
  DT_mn = DT[, list(mean(value)), by = c("run", "sim_yr", "wym")]
  DT_mn = merge.data.table(DT_mn, scenarios, by = "run")
  DT_mn$yr_mn = (DT_mn$wym + (DT_mn$sim_yr - min(DT_mn$sim_yr)) * 12)
  return(DT_mn)
}


mn_data = list()
p = list()

vars = c("plantc", "lai", "height", "psn")

for (i in seq_along(vars)) {
  path = paste0("output/spinup/R_obj/", vars[i], ".rdata")
  DT_all = read_robj(path)
  mn_data[[i]] = mean2mn(DT_all, spin_scenarios)
  
  p[[i]] =  ggplot(mn_data[[i]]) +
    aes(x = yr_mn, y = V1, color = veg) +
    geom_line(size = 1L) +
    labs(x = "Month",
         y = vars[i],
         title = paste0("Monthly", vars[i])) +
    theme_light()
  
}

grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], ncol = 2)





