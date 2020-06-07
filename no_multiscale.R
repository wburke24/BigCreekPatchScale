# no_multiscale

setwdhere()

library(RHESSysIOinR)


# -------------------- Multiscale dev stuff --------------------
comp_ID = Sys.info()["nodename"]
compile = FALSE
if (compile) {
  comp_out = compile_rhessys(location = "../../../rhessys/RHESSys/rhessys/", delete_objs = FALSE, 
                             destination = file.path("bin",comp_ID), CFLAGS = "-w -O2")
}



# s facing conifer
worldfile = "worldfiles/no_multi/BC_conifer_S_1p.world"

# redefs
# thin 0.4
# build_redefine(worldfile, paste0(worldfile, ".Y1944M10D1H10"),  std_thin = "0.6", veg_parm_ID = "50")

# prescribed fire
# litter_vars = c("litter_cs.litr1c", "litter_ns.litr1n",	"litter_cs.litr2c", "litter_cs.litr3c", "litter_cs.litr4c", "cs.cwdc", "ns.cwdn")
# build_redefine(worldfile, paste0(worldfile, ".Y1944M11D1H10"), vars = litter_vars, values = 0)


# -------------------- Computer ID --------------------
comp_ID = Sys.info()["nodename"]

# load scenarios
# scenarios = NULL
# load("scenarios.rdata")

# -------------------- Project/Run Name --------------------
name = "BC_patch_nomult"
n_sim = 1

# -------------------- Climate Data --------------------
clim_data = "clim/upGGmod"
start_end = c("1943 10 01 01", "1963 10 01 01")

# -------------------- Parameter Method --------------------
parameter_method = "exact_values"

# -------------------- Input RHESSys --------------------
input_rhessys <- list()
#input_rhessys$rhessys_version <- file.path("bin",comp_ID,"rhessys7.1.1")
input_rhessys$rhessys_version <- file.path("bin",comp_ID,"rhessys7.1")
input_rhessys$tec_file <- "tecfiles/BC.tec"
input_rhessys$world_file <- "worldfiles/no_multi/BC_conifer_S_1p.world"
input_rhessys$world_hdr_prefix <- "BC"
input_rhessys$flow_file <- "flowtables/BC_1p.flow"
input_rhessys$start_date <- start_end[1]
input_rhessys$end_date <- start_end[2]
input_rhessys$output_folder <- "output/"
input_rhessys$output_filename <- "BC_1p"
input_rhessys$command_options <- c("-g -b -p -c") #-v -6  

# -------------------- Input Headers --------------------
input_hdr_list <- list()
input_hdr_list$basin_def <- c("defs/basin_p301.def")
input_hdr_list$hillslope_def <- c("defs/hill_p301.def")
input_hdr_list$zone_def <- c("defs/zone_p301.def")
input_hdr_list$soil_def <- c("defs/soil_forestshrub.def")
#input_hdr_list$soil_def <- c("defs/soil_sandyloam.def")
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
input_standard_par_list <- list(
  m = c(2),
  k = c(2),
  m_v = c(2),
  k_v = c(2),
  pa = c(1.15),
  po = c(0.766),
  gw1 = c(0.24),
  gw2 = c(0.2)
)

# -------------------- Make Climate Basestation --------------------
#input_clim_base_list <- NULL
input_clim_base_list = clim_auto(base_station_id = 101,
                                 x_coordinate = 100.0,
                                 y_coordinate = 100.0,
                                 z_coordinate = 1748,
                                 effective_lai = 3.5,
                                 screen_height = 2,
                                 daily_prefix = "clim/upGGmod")

# -------------------- Make Dated Sequence --------------------
input_dated_seq_list <- NULL

# -------------------- Make Tec File --------------------
# tec handled inside iteration loop
input_tec_data = input_tec(NULL, start_end = start_end)


# -------------------- Output --------------------
output_method = "r"
# old vars at basin level
# output_variables <- data.frame(out_file = character(), variable = character(), stringsAsFactors = FALSE)
# output_variables[1,] <- data.frame("bd", "streamflow", stringsAsFactors = FALSE)
# output_variables[2,] <- data.frame("bd", "plantc", stringsAsFactors = FALSE)
# output_variables[3,] <- data.frame("bd", "psn", stringsAsFactors = FALSE)
# output_variables[4,] <- data.frame("bd", "lai", stringsAsFactors = FALSE)
# output_variables[5,] <- data.frame("bd", "litter_capacity", stringsAsFactors = FALSE)
# output_variables[6,] <- data.frame("bd", "litter_store", stringsAsFactors = FALSE)

# ----- veg, soil and sharing parameters -----
input_def_list = soilslist$medium

def_df <- data.frame(matrix(unlist(input_def_list), nrow = length(input_def_list), byrow = T), stringsAsFactors = FALSE)

# this one has 2 treatments - thin from below 
# and prescribed fire
# tdate1 = "1944 10 01 01"
# tec_treat = tec_repeat(start = tdate1, end = start_end[2], interval = NA, unit = NA, event_name = "redefine_world_thin_remain")
# tdate2 = "1944 11 01 01"
# tec_treat2 = tec_repeat(start = tdate2, end = start_end[2], interval = NA, unit = NA, event_name = "redefine_world")
# 
# input_tec_data = rbind(input_tec_data, tec_treat)
# input_tec_data = rbind(input_tec_data, tec_treat, tec_treat2)
# input_tec_data = input_tec_data[with(input_tec_data, order(year, month, day, hour)), ]
# 

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
# data_out = select_output_variables_R(output_variables = output_variables,
#                                      output_folder = input_rhessys$output_folder,
#                                      output_filename = input_rhessys$output_filename,
#                                      run = 1,
#                                      max_run =  1,
#                                      return_data = TRUE)



BC_1p = readin_rhessys_output(pre = "output/BC_1p")

BC_1p_cdg = BC_1p$cdg
BC_1p_cd = BC_1p$cd
BC_1p_pdg = BC_1p$pdg
BC_1p_pd = BC_1p$pd

library(ggplot2)
library(colorspace)

ggplot(BC_1p_pdg) +
  aes(x = date) +
  geom_line(aes(y = litr1c), size = 1L, colour = rainbow_hcl(3)[1]) +
  geom_line(aes(y = litr2c), size = 1L, colour = rainbow_hcl(3)[2]) +
  geom_line(aes(y = litr3c), size = 1L, colour = rainbow_hcl(3)[3]) +
  theme_classic()

BC_1p_cdg$stratumID = as.factor(BC_1p_cdg$stratumID)
BC_1p_cd$stratumID = as.factor(BC_1p_cd$stratumID)


ggplot(BC_1p_cdg) +
 aes(x = date, y = leafc, colour = stratumID) +
 geom_line(size = 1L) +
 scale_color_hue() +
 theme_classic()

ggplot(BC_1p_cd) +
 aes(x = date, y = height, colour = stratumID) +
 geom_line(size = 1L) +
 scale_color_hue() +
 theme_classic()

ggplot(BC_1p_cdg) +
 aes(x = date, y = cwdc, colour = stratumID) +
 geom_line(size = 1L) +
 scale_color_hue() +
 theme_classic()



