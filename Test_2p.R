setwdhere()

library(RHESSysIOinR)

# -------------------- Multiscale dev stuff --------------------
comp_ID = Sys.info()["nodename"]
compile = FALSE
if (compile) {
  comp_out = compile_rhessys(location = "../../../rhessys/rhessys-develop/rhessys/",
                             destination = file.path("bin",comp_ID))
}

# -------------------- Project/Run Name --------------------
# name = "BC_2p"

# -------------------- Climate Data --------------------
# clim_data = "clim/upGGmod"
# start_end = c("1943 10 01 01", "1948 10 01 01")

# -------------------- Input RHESSys --------------------
input_rhessys = IOin_rhessys_input(
  version = file.path("bin", comp_ID, "rhessys7.2"),
  tec_file = "tecfiles/BC.tec",
  world_file = "worldfiles/no_multi/BC_conifer_S_2p.world",
  world_hdr_prefix = "BC",
  flowtable = "flowtables/BC_2p.flow",
  start = "1943 10 01 01",
  end = "1948 10 01 01",
  output_folder = "output/",
  output_prefix = "BC_2p",
  commandline_options = c("-vmort_off -g -b -p -c -msr -v -6")
)

# -------------------- Input Headers --------------------
input_hdr_list = IOin_hdr(
  basin = "defs/basin_p301.def",
  hillslope = "defs/hill_p301.def",
  zone = "defs/zone_p301.def",
  soil = "defs/soil_forestshrub.def",
  landuse = "defs/lu_p301.def",
  stratum = c("defs/veg_p301_shrub.def", "defs/veg_p301_conifer.def"),
  basestations = "clim/upGGmod.base"
)

# --------------------  Def File Parameters --------------------
#input_def_list = NULL

# -------------------- Soils --------------------
soilslist = NULL
load("soilslist.rdata")

input_def_pars = IOin_def_pars_simple(
  list("defs/veg_p301_conifer.def", "epc.allocation_flag","dickenson"),
  list("defs/veg_p301_conifer.def", "epc.alloc_frootc_leafc", 1),
  list("defs/veg_p301_conifer.def", "epc.alloc_stemc_leafc", 0.6),
  list("defs/veg_p301_conifer.def", "epc.netpabs_shade", 0.2),
  list("defs/veg_p301_conifer.def", "epc.netpabs_sunlit", 0.2),
  list("defs/veg_p301_conifer.def", "epc.root_growth_direction", 0.9),
  list("defs/veg_p301_conifer.def", "epc.proj_sla", 8),
  list("defs/veg_p301_conifer.def", "epc.gl_smax", 0.005),
  list("defs/veg_p301_conifer.def", "epc.resprout_leaf_carbon", 0.01),
  list("defs/veg_p301_conifer.def", "epc.min_percent_leafg", 0.01),
  list("defs/veg_p301_conifer.def", "epc.min_leaf_carbon", 0.001),
  list("defs/veg_p301_conifer.def", "epc.tcoef", 0.6),
  list("defs/veg_p301_conifer.def", "epc.topt", 15),
  list("defs/veg_p301_conifer.def", "epc.branch_turnover", 0.005),
  list("defs/veg_p301_conifer.def", "epc.max_daily_mortality", 0.01),
  list("defs/veg_p301_conifer.def", "epc.min_daily_mortality", 0.01),
  list("defs/veg_p301_conifer.def", "epc.froot_turnover", 0.5),
  list("defs/veg_p301_conifer.def", "epc.leaf_turnover", 0.5),
  list("defs/veg_p301_conifer.def", "epc.gr_perc", 0.18),
  list("defs/veg_p301_conifer.def", "mrc.per_N", 0.21),
  list("defs/veg_p301_conifer.def", "epc.ext_coef", 0.7),
  list("defs/veg_p301_conifer.def", "epc.storage_transfer_prop", 0.7),
  list("defs/veg_p301_conifer.def", "epc.root_distrib_parm", 7),
  list("defs/soil_forestshrub.def", "soil_depth", 1),
  list("defs/soil_forestshrub.def", "pore_size_index", 0.15),
  list("defs/soil_forestshrub.def", "sat_to_gw_coeff", 0.12),
  list("defs/soil_forestshrub.def", "psi_air_entry", 0.15)
)

input_def_pars_sets = IOin_def_pars_simple(input_def_pars, n = 10)



# -------------------- Standard (soil/subsurface) Parameters--------------------
input_std_pars <- list(m = c(2), k = c(2), m_v = c(2), k_v = c(2), pa = c(1.15),
  po = c(0.766), gw1 = c(0.24), gw2 = c(0.2)
)

std_pars = IOin_std_pars(
  m = 2,
  k = 2,
  m_v = 2,
  k_v = 2,
  pa = 1.15,
  po = 0.766,
  gw1 = 0.24,
  gw2 = 0.2
)

# -------------------- Make Tec File --------------------
input_tec_data = IOin_tec_std(start_end = start_end)

# ----- veg, soil and sharing parameters -----
#input_def_list = soilslist$medium

# -------------------- End Run Setup--------------------

# -------------------- Run RHESSys --------------------

run_rhessys_core(
  input_rhessys = input_rhessys,
  hdr_files = input_hdr_list,
  def_pars = input_def_pars,
  std_pars = input_std_pars,
  tec_data = input_tec_data
)

rhout = readin_rhessys_output("output/BC_2p")

wb = watbal_basin(bd = rhout$bd)

summary(wb$watbal)


