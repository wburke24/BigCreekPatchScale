# soils_options

input_hdr_list <- list()
input_hdr_list$basin_def <- c("defs/basin_p301.def")
input_hdr_list$hillslope_def <- c("defs/hill_p301.def")
input_hdr_list$zone_def <- c("defs/zone_p301.def")
input_hdr_list$soil_def <- c("defs/soil_forestshrub.def")
#input_hdr_list$soil_def <- c("defs/soil_sandyloam.def")
input_hdr_list$landuse_def <- c("defs/lu_p301.def")
input_hdr_list$patch_def <- c("defs/patch_p301.def")
input_hdr_list$stratum_def <- c("defs/veg_p301_shrub.def", "defs/veg_p301_conifer.def")
input_hdr_list$base_stations <- c("clim/Grove_lowprov_clim.base")

# --- Shallow soils ---
shallow <- list(
  list(input_hdr_list$stratum_def[2], "epc.allocation_flag", c("dickenson")), #"../defs/4ponderosa.pine.comp"
  list(input_hdr_list$stratum_def[2], "epc.alloc_frootc_leafc", 1.5 ),
  list(input_hdr_list$stratum_def[2], "epc.alloc_stemc_leafc", 0.6),
  list(input_hdr_list$stratum_def[2], "epc.netpabs_shade" , 0.2),
  list(input_hdr_list$stratum_def[2], "epc.netpabs_sunlit" , 0.2),
  list(input_hdr_list$stratum_def[2], "epc.psi_close", -9.0),
  list(input_hdr_list$stratum_def[2], "epc.root_growth_direction",0.85),
  list(input_hdr_list$stratum_def[2], "epc.proj_sla",7),
  list(input_hdr_list$stratum_def[2], "epc.gl_smax", 0.004),
  list(input_hdr_list$stratum_def[2], "epc.min_percent_leafg", 0.01),
  list(input_hdr_list$stratum_def[2], "epc.resprout_leaf_carbon", 0.01),
  list(input_hdr_list$stratum_def[2], "epc.min_leaf_carbon", 0.01),
  list(input_hdr_list$stratum_def[2], "epc.tcoef", 0.6),
  list(input_hdr_list$stratum_def[2], "epc.topt", 15.0),
  list(input_hdr_list$stratum_def[2], "epc.branch_turnover", 0.006),
  list(input_hdr_list$stratum_def[2], "epc.max_daily_mortality", 0.01),
  list(input_hdr_list$stratum_def[2], "epc.min_daily_mortality", 0.01),
  list(input_hdr_list$stratum_def[2], "epc.froot_turnover", 0.3),
  list(input_hdr_list$stratum_def[2], "epc.leaf_turnover", 0.4),
  list(input_hdr_list$stratum_def[2], "epc.gr_perc", 0.18),
  list(input_hdr_list$stratum_def[2], "mrc.per_N", 0.21),
  list(input_hdr_list$stratum_def[2], "epc.ext_coef", 0.7),
  list(input_hdr_list$stratum_def[2], "epc.storage_transfer_prop", 0.5),
  list(input_hdr_list$stratum_def[2], "epc.root_distrib_parm", 10.000000), #../defs/4ponderosa.pine.tmp)
  list(input_hdr_list$soil_def[1], "soil_depth", 0.25), #../defs/patch_fp301.def > tmp2)
  list(input_hdr_list$soil_def[1], "pore_size_index", 0.15),
  list(input_hdr_list$soil_def[1], "sat_to_gw_coeff", 0.2),
  list(input_hdr_list$soil_def[1], "psi_air_entry", 0.16) #../defs/patch_fp301.tmp)
)
# --- Medium soils ---
medium <- list(
  list(input_hdr_list$stratum_def[2], "epc.allocation_flag", "dickenson"), # ../defs/4ponderosa.pine.comp
  list(input_hdr_list$stratum_def[2], "epc.alloc_frootc_leafc",1.0),
  list(input_hdr_list$stratum_def[2], "epc.alloc_stemc_leafc", 0.6),
  list(input_hdr_list$stratum_def[2], "epc.netpabs_shade", 0.2),
  list(input_hdr_list$stratum_def[2], "epc.netpabs_sunlit", 0.2),
  list(input_hdr_list$stratum_def[2], "epc.root_growth_direction", 0.9),
  list(input_hdr_list$stratum_def[2], "epc.proj_sla", 8),
  list(input_hdr_list$stratum_def[2], "epc.gl_smax", 0.005),
  list(input_hdr_list$stratum_def[2], "epc.resprout_leaf_carbon", 0.01),
  list(input_hdr_list$stratum_def[2], "epc.min_percent_leafg", 0.01),
  list(input_hdr_list$stratum_def[2], "epc.min_leaf_carbon", 0.001),
  list(input_hdr_list$stratum_def[2], "epc.tcoef", 0.6),
  list(input_hdr_list$stratum_def[2], "epc.topt", 15.0),
  list(input_hdr_list$stratum_def[2], "epc.branch_turnover", 0.005),
  list(input_hdr_list$stratum_def[2], "epc.max_daily_mortality", 0.01),
  list(input_hdr_list$stratum_def[2], "epc.min_daily_mortality", 0.01),
  list(input_hdr_list$stratum_def[2], "epc.froot_turnover", 0.5),
  list(input_hdr_list$stratum_def[2], "epc.leaf_turnover", 0.5),
  list(input_hdr_list$stratum_def[2], "epc.gr_perc", 0.18),
  list(input_hdr_list$stratum_def[2], "mrc.per_N", 0.21),
  list(input_hdr_list$stratum_def[2], "epc.ext_coef", 0.7),
  list(input_hdr_list$stratum_def[2], "epc.storage_transfer_prop", 0.7),
  list(input_hdr_list$stratum_def[2], "epc.root_distrib_parm", 7.000000), # ../defs/4ponderosa.pine.tmp
  list(input_hdr_list$soil_def[1], "soil_depth", 1.0), # ../defs/patch_fp301.def
  list(input_hdr_list$soil_def[1], "pore_size_index", 0.15),
  list(input_hdr_list$soil_def[1], "sat_to_gw_coeff", 0.12),
  list(input_hdr_list$soil_def[1], "psi_air_entry", 0.15) # ../defs/patch_fp301.tmp
)
# --- Deep soils ---
deep <- list(
  list(input_hdr_list$stratum_def[2], "epc.allocation_flag", "dickenson"), # ../defs/4ponderosa.pine.comp 
  list(input_hdr_list$stratum_def[2], "epc.alloc_frootc_leafc", 1.0),
  list(input_hdr_list$stratum_def[2], "epc.alloc_stemc_leafc", 0.6),
  list(input_hdr_list$stratum_def[2], "epc.netpabs_shade", 0.35),
  list(input_hdr_list$stratum_def[2], "epc.netpabs_sunlit", 0.35),
  list(input_hdr_list$stratum_def[2], "epc.root_growth_direction", 0.95),
  list(input_hdr_list$stratum_def[2], "epc.proj_sla", 8),
  list(input_hdr_list$stratum_def[2], "epc.gl_smax", 0.005),
  list(input_hdr_list$stratum_def[2], "epc.resprout_leaf_carbon", 0.05),
  list(input_hdr_list$stratum_def[2], "epc.min_percent_leafg", 0.01),
  list(input_hdr_list$stratum_def[2], "epc.min_leaf_carbon", 0.001),
  list(input_hdr_list$stratum_def[2], "epc.tcoef", 0.6),
  list(input_hdr_list$stratum_def[2], "epc.topt", 15.0),
  list(input_hdr_list$stratum_def[2], "epc.branch_turnover", 0.005),
  list(input_hdr_list$stratum_def[2], "epc.max_daily_mortality", 0.01),
  list(input_hdr_list$stratum_def[2], "epc.min_daily_mortality", 0.01),
  list(input_hdr_list$stratum_def[2], "epc.froot_turnover", 0.5),
  list(input_hdr_list$stratum_def[2], "epc.leaf_turnover", 0.5),
  list(input_hdr_list$stratum_def[2], "epc.gr_perc", 0.18),
  list(input_hdr_list$stratum_def[2], "mrc.per_N", 0.21),
  list(input_hdr_list$stratum_def[2], "epc.ext_coef", 0.7),
  list(input_hdr_list$stratum_def[2], "epc.storage_transfer_prop", 0.7),
  list(input_hdr_list$stratum_def[2], "epc.root_distrib_parm", 2.000000), # ../defs/4ponderosa.pine.tmp 
  list(input_hdr_list$soil_def[1], "soil_depth", 6.0), # ../defs/patch_fp301.def
  list(input_hdr_list$soil_def[1], "pore_size_index", 0.2),
  list(input_hdr_list$soil_def[1], "sat_to_gw_coeff", 0.08),
  list(input_hdr_list$soil_def[1], "psi_air_entry", 0.218) # ../defs/patch_fp301.tmp
)

soilslist = list(shallow, medium, deep)
names(soilslist) = c("shallow", "medium", "deep")

save(soilslist,file = "soilslist.rdata")
