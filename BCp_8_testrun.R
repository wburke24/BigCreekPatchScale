# -------------------- Start --------------------
setwdhere()

repull = FALSE
if (repull) {
  devtools::install_github("RHESSys/RHESSysIOinR")
  devtools::install_github("RHESSys/RHESSysPreprocessing")
  devtools::install("../../../rhessys/RHESSysIOinR")
}

library(RHESSysIOinR)
library(RHESSysPreprocessing)
library(data.table)
library(tictoc)


options(stringsAsFactors = FALSE) # just in case

# -------------------- Computer ID --------------------
comp_ID = Sys.info()["nodename"]

# -------------------- Multiscale dev stuff --------------------
compile = FALSE
if (compile) {
  comp_out = compile_rhessys(location = "../../../rhessys/rhessys-multiscale", delete_objs = TRUE, 
                             destination = file.path("bin",comp_ID), CFLAGS = "-w -O2")
}

# --- load scenarios ---
scenarios = NULL
load("scenarios.rdata")

# select some limited scenarios
# - all for conifer over shrub
# - thin from below at .4, selection thin w removal at .4, no treatment
#  - 5 and 20 yr interval
# - deep soils, 0.5 sharing, S aspect
# - all clims, 

scenarios = data.table(scenarios)
# scenarios = scenarios[veg == "conifer" & aspect == "N" & sharing == 0.5 & 
#                         ((interval == "none" & treatments == "no treatment") | 
#                            (interval %in% c("5", "20") & treatments %in% c("over_0.4 carbon_removed", "under_0.4 1mo_rm100litter") ) ) , ]
# 
# scenarios = scenarios[veg == "conifer" & aspect == "N" & sharing == 0.5 & interval == "20" & treatments == "over_0.4 carbon_removed" & 
#                         climate == "clim/upGGmod" & start_period == "variable" & soils == "deep"  , ]


# scenarios = scenarios[veg == "shrub" & aspect == "N" & sharing == 0.5 & climate == "clim/upGGmod" & start_period == "variable" &
#                         ((interval == "none" & treatments == "no treatment") | 
#                            (interval %in% c("5", "20") & treatments %in% c("under_0.4 1mo_rm100litter") ) ) , ]

scenarios = scenarios[veg == "shrub" & aspect == "N" & sharing == 0.5 & interval == "20" & treatments == "under_0.4 1mo_rm100litter" &
                        climate == "clim/upGGmod" & start_period == "variable" & soils == "deep"  , ]


#(interval == "none" & treatments == "no treatment")
scenarios$end_date = "2010 10 01 01"



# -------------------- Project/Run Name --------------------
name = "BCpatch_noshade"
n_sim = 1

# -------------------- Parameter Method --------------------
parameter_method = "exact_values"

# -------------------- Input RHESSys --------------------
input_rhessys <- list()
input_rhessys$rhessys_version <- file.path("bin",comp_ID,"rhessys7.1.1")
input_rhessys$tec_file <- "tecfiles/BC.tec"
input_rhessys$world_file <- NULL
input_rhessys$world_hdr_prefix <- "BC"
input_rhessys$flow_file <- NULL
input_rhessys$start_date <- NULL
input_rhessys$end_date <- NULL
input_rhessys$output_folder <- "output/"
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
input_hdr_list$base_stations <- c("clim/upGGmod.base")

# -------------------- Parameters --------------------
input_preexisting_table <- NULL

# --------------------  Def File Parameters --------------------
#input_def_list = NULL

# -------------------- Soils --------------------
soilslist = NULL
load("soilslist.rdata")

# -------------------- Standard (soil/subsurface) Parameters--------------------
input_standard_par_list <- list(m = c(2), k = c(2), m_v = c(2), k_v = c(2), pa = c(1.15), po = c(0.766), gw1 = c(0.24), gw2 = c(0.2))

# -------------------- Make Climate Basestation --------------------
#input_clim_base_list <- NULL
input_clim_base_list = clim_auto(base_station_id = 101, x_coordinate = 100.0, y_coordinate = 100.0, z_coordinate = 1748, 
                                 effective_lai = 3.5, screen_height = 2, daily_prefix = "clim/upGGmod")

# -------------------- Make Dated Sequence --------------------
input_dated_seq_list <- NULL

# -------------------- Make Tec File --------------------
# tec handled inside iteration loop
# input_tec_data = input_tec(NULL, start_end = start_end)

# -------------------- Output --------------------
output_method = "r"

# vars at basin level
# output_variables <- data.frame(out_file = character(), variable = character(), stringsAsFactors = FALSE)
# output_variables[1,] <- data.frame("bd", "streamflow", stringsAsFactors = FALSE)
# output_variables[2,] <- data.frame("bd", "plantc", stringsAsFactors = FALSE)
# output_variables[3,] <- data.frame("bd", "psn", stringsAsFactors = FALSE)
# output_variables[4,] <- data.frame("bd", "lai", stringsAsFactors = FALSE)
# output_variables[5,] <- data.frame("bd", "litter_capacity", stringsAsFactors = FALSE)
# output_variables[6,] <- data.frame("bd", "litter_store", stringsAsFactors = FALSE)

# patch/canopy level vars
output_variables <- data.frame(out_file = character(), variable = character(), stringsAsFactors = FALSE)
output_variables[1,] <- data.frame("pd", "Qout")
output_variables[2,] <- data.frame("pd", "litter.S") # fuel/litter moisture
output_variables[3,] <- data.frame("pd", "psn")
output_variables[4,] <- data.frame("pdg", "lai")
output_variables[5,] <- data.frame("pdg", "plantc")
output_variables[6,] <- data.frame("cd", "height")
output_variables[7,] <- data.frame("cdg", "leafc") # add fuel veg/density from maureens: cover frac * cs.leafc
# fuel litter: litter_cs.litr1c +	litter_cs.litr2c +	litter_cs.litr3c +	litter_cs.litr4c + cs.dead_leafc
output_variables[8,] <- data.frame("pdg", "litrc_sum") # sum of litr C
output_variables[9,] <- data.frame("pd", "rz_storage") # soil moisture: rootzone.S
output_variables[10,] <- data.frame("pd", "fire_et")
output_variables[11,] <- data.frame("pd", "pet")

# other fire vars - ingoring for now: wind, wind_direction, relative_humidity, z, temp
#output_variables$id_extract = NULL

# ========================= Iterate Scenarios =========================



# ||||| -- RUN FROM HERE -- |||||

tic()


veg = unique(scenarios$veg)

for (v in seq_along(veg)) {
  
  runs = which(scenarios$veg == veg[v])
  
  input_rhessys$output_filename = veg[v]
  
  if (!dir.exists("output/allsim")) {
    dir.create("output/allsim")
  }
  
  # run count for output
  runct = 0
  
  for (i in runs) {
    
    runct = runct + 1
    cat("\n ========== Scenario ",runct," ========== \n\n")
    
    # ----- input_rhessys vars -----
    input_rhessys$start_date <- scenarios$start_date[i]
    input_rhessys$end_date <- scenarios$end_date[i]
    input_rhessys$world_file = scenarios$worldfile[i]
    input_rhessys$flow_file = paste0("flowtables/BC_", scenarios$veg[i], ".flow")
    
    # ----- veg, soil and sharing parameters -----
    input_def_list = c(list(list(input_hdr_list$landuse_def[1], "sh_l", scenarios$sharing[i]),
                            list(input_hdr_list$landuse_def[1], "shading_flag", 0),
                            list(input_hdr_list$stratum_def[1], "epc.height_to_stem_coef", 2.0),
                            list(input_hdr_list$stratum_def[2], "epc.height_to_stem_coef", 6.0)),
                       soilslist[[scenarios$soils[i]]] )
    
    # ----- climate/CO2 -----
    input_clim_base_list[[1]]$daily$c1[1] = scenarios$climate[i]
    if (scenarios$climate[i] == "clim/upGGmod_2C") {
      input_def_list = c(input_def_list, list(list(input_hdr_list$zone_def[1], "atm_CO2", 450)))
    }
    # ----- treatments - tec events -----
    input_tec_data = input_tec(start = scenarios$start_date[i], end = scenarios$end_date[i], output_state = FALSE)
    
    # thinning type - NA is no thinning (could still have prescribed fire)
    if (!is.na(scenarios$treat_type[i])) {
      redef_world = paste0("BC_", scenarios$veg[i], "_", scenarios$treat_type[i], "_", scenarios$treat_intensity[i], ".world")
      redef_path = file.path("worldfiles/redefine", redef_world)
      tec_treat = tec_repeat(scenarios$start_date[i], scenarios$end_date[i], scenarios$interval[i], "year", 
                             paste0("redefine_world_", scenarios$tec_type[i]), scenarios$worldfile[i], redef_path, overwrite = TRUE)
      input_tec_data = rbind(input_tec_data, tec_treat)
    }
    if (scenarios$presc_fire[i]) {
      redef_world = paste0("BC_", scenarios$veg[i], "_litter_0.world")
      redef_path = file.path("worldfiles/redefine", redef_world)
      start_1mo = paste0(gsub("-"," ", seq(as.Date(scenarios$start_date[i], "%Y %m %d %H"), by = "month", length = 2)[2]), " 10")
      tec_presc = tec_repeat(start_1mo, scenarios$end_date[i], scenarios$interval[i], "year", "redefine_world", 
                             scenarios$worldfile[i], redef_path, overwrite = TRUE)
      input_tec_data = rbind(input_tec_data, tec_presc)
    }
    input_tec_data = input_tec_data[with(input_tec_data, order(year, month, day, hour)), ]
    
    # -------------------- End Run Setup--------------------

    input_tec_data = input_tec_data[1:3,]
    
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
    
    # ----- clean up redefines -----
    shh = file.remove(file.path(dirname(scenarios$worldfile[i]), 
                                list.files(path = dirname(scenarios$worldfile[i]), include.dirs = FALSE, pattern = "H10$|H1$" )))
    
    # -------------------- Get Output --------------------
    
    #substr(scenarios$worldfile[i], 12, nchar(scenarios$worldfile[i]) - 6)
    selout = F
    if (selout) {
      data_out = select_output_variables_R(output_variables = output_variables,
                                           output_folder = input_rhessys$output_folder,
                                           output_filename = input_rhessys$output_filename,
                                           run = runct,
                                           max_run =  length(runs),
                                           return_data = FALSE)
    }
    
    cat("\n ========== End Scenario ",runct," ========== \n\n")
    
  } # end runs loop
  
  file.rename("output/allsim", paste0("output/", veg[v]))
  
} # end veg loop

toc()

# ==================== end scenario iterations ====================



veg = unique(scenarios$veg)
folders = paste0("output/", veg)
# get from run script 
vars = output_variables$variable

for (i in vars) {
  
  DT_all = NULL
  
  for (j in seq_along(folders)) {
    
    # get the data
    read_in = fread(file.path(folders[j], i)) 
    
    # set colnames to just the nums
    colnames(read_in)[startsWith(colnames(read_in), "run")] = which(scenarios$veg == veg[j])
    
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
    dt_melt[, run := as.integer(run), ]
    #dt_melt$run = as.integer(dt_melt$run)
    #colnames(dt_melt)[colnames(dt_melt) == "value"] = i
    
    if (j == 1) {
      DT_all = dt_melt
    } else {
      DT_all = rbind(DT_all, dt_melt)
    }
    
  }
  
  save(DT_all, file = paste0("output/R_obj_noshade/", i, ".rdata"))
  
  rm(read_in, dt, dt_agg, dt_melt)
  gc()
  
}




