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
  comp_out = compile_rhessys(location = "../../../rhessys/rhessys-develop", delete_objs = TRUE, 
                             destination = file.path("bin",comp_ID), CFLAGS = "-w -O2")
}

# --- load scenarios ---
scenarios = NULL
load("scenarios.rdata")

scenarios = as.data.table(scenarios)
scenarios = scenarios[c(3467L), ]

# -------------------- Project/Run Name --------------------
name = "BC_patch"
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
                                 effective_lai = 3.5, screen_height = 2, daily_prefix = "clim/upGG_39_15")

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
output_variables[1,] = c("pd", "Qout")
output_variables[2,] <- c("pd", "litter.S") # fuel/litter moisture
output_variables[3,] <- c("pd", "psn")
output_variables[4,] <- c("pdg", "lai")
output_variables[5,] <- c("pdg", "plantc")
output_variables[6,] <- c("cd", "height")
output_variables[7,] <- c("cdg", "leafc") # add fuel veg/density from maureens: cover frac * cs.leafc
# fuel litter: litter_cs.litr1c +	litter_cs.litr2c +	litter_cs.litr3c +	litter_cs.litr4c + cs.dead_leafc
output_variables[8,] <- c("pdg", "litrc_sum") # sum of litr C
output_variables[9,] <- c("pd", "rz_storage") # soil moisture: rootzone.S
output_variables[10,] <- c("pd", "fire_et")
output_variables[11,] <- c("pd", "pet")

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
     
    # SET START TO -3 YEARS
    sim_start = paste0(as.numeric(substr(scenarios$start_date[i],0,4)) - 3, substr(scenarios$start_date[i], 5, nchar(scenarios$start_date[i]) ) )
    input_rhessys$start_date <- sim_start
    input_rhessys$end_date <- scenarios$end_date[i]
    input_rhessys$world_file = scenarios$worldfile[i]
    input_rhessys$flow_file = paste0("flowtables/BC_", scenarios$veg[i], ".flow")
    
    # ----- veg, soil and sharing parameters ----- stratum 1 = shrub, stratum 2 = conifer
    input_def_list = c(list(list(input_hdr_list$landuse_def[1], "sh_l", scenarios$sharing[i]),
                            list(input_hdr_list$stratum_def[1], "epc.height_to_stem_coef", 2.0),
                            list(input_hdr_list$stratum_def[1], "epc.resprout_leaf_carbon", 0.01),                       
                            list(input_hdr_list$stratum_def[2], "epc.height_to_stem_coef", 6.0),
                            list(input_hdr_list$stratum_def[2], "epc.resprout_leaf_carbon", 0.05),
                            list(input_hdr_list$stratum_def[2], "epc.zero_turnover_sprouts", 1)),
                       soilslist[[scenarios$soils[i]]] )
    
    # ----- climate/CO2 -----
    input_clim_base_list[[1]]$daily$c1[1] = scenarios$climate[i]
    if (scenarios$climate[i] == "clim/upGG_39_15_2C") {
      input_def_list = c(input_def_list, list(list(input_hdr_list$zone_def[1], "atm_CO2", 450)))
    }
    # ----- treatments - tec events -----
    input_tec_data = input_tec(start = scenarios$start_date[i], end = scenarios$end_date[i], output_state = FALSE)
    
    # thinning type - NA is no thinning (could still have prescribed fire)
    if (!is.na(scenarios$treat_type[i])) {
      redef_world = paste0("BC_", scenarios$veg[i], "_", scenarios$treat_type[i], "_", scenarios$treat_intensity[i], ".world")
      redef_path = file.path("worldfiles/redefine", redef_world)
      
      #CHANGE DATE HERE
      treat_date = paste0(as.numeric(substr(scenarios$start_date[i],0,4)) + 1, " 04 01 01" )
      tec_treat = tec_repeat(treat_date, scenarios$end_date[i], scenarios$interval[i], "year", 
                             paste0("redefine_world_", scenarios$tec_type[i]), scenarios$worldfile[i], redef_path, overwrite = TRUE)
      input_tec_data = rbind(input_tec_data, tec_treat)
    }
    if (scenarios$presc_fire[i]) {
      redef_world = paste0("BC_", scenarios$veg[i], "_litter_0.world")
      redef_path = file.path("worldfiles/redefine", redef_world)
      treat_date2 = paste0(as.numeric(substr(scenarios$start_date[i],0,4)) + 1, " 05 01 01" )
      #start_1mo = paste0(gsub("-"," ", seq(as.Date(scenarios$start_date[i], "%Y %m %d %H"), by = "month", length = 2)[2]), " 10")
      tec_presc = tec_repeat(treat_date2, scenarios$end_date[i], scenarios$interval[i], "year", "redefine_world", 
                             scenarios$worldfile[i], redef_path, overwrite = TRUE)
      input_tec_data = rbind(input_tec_data, tec_presc)
    }
    input_tec_data = input_tec_data[with(input_tec_data, order(year, month, day, hour)), ]
    
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
    
    # ----- clean up redefines -----
    shh = file.remove(file.path(dirname(scenarios$worldfile[i]), 
                                list.files(path = dirname(scenarios$worldfile[i]), include.dirs = FALSE, pattern = "H10$|H1$" )))
    
    # -------------------- Get Output --------------------
    
    #substr(scenarios$worldfile[i], 12, nchar(scenarios$worldfile[i]) - 6)
    
    data_out = select_output_variables_R(output_variables = output_variables,
                                         output_folder = input_rhessys$output_folder,
                                         output_filename = input_rhessys$output_filename,
                                         run = runct,
                                         max_run =  length(runs),
                                         return_data = FALSE)
    
    cat("\n ========== End Scenario ",runct," ========== \n\n")
    
  } # end runs loop
  
  if (file.exists(paste0("output/", veg[v]))) {
    file.rename(paste0("output/", veg[v]), paste0("output/", veg[v],"_old"))
  }
  file.rename("output/allsim", paste0("output/", veg[v]))
  
} # end veg loop

toc()

# ==================== end scenario iterations ====================

cleanup_rhessys(dir = "output/")



