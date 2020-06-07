# BCp_init_worldfiles

setwdhere()

#devtools::install_github("RHESSys/RHESSysPreprocessing")
#devtools::install("../rhessys/RHESSysPreprocessing")
library(RHESSysPreprocessing)
library(data.table)

# ----- get the values from awks -----
get_vars = function(path) {
  inlines = readLines(path)
  
  var_match = regmatches(inlines, gregexpr("\\$2 == \\\".*\\\")", inlines))
  var_chr = unlist(var_match[var_match != "character(0)"])
  vars = gsub(".{2}$","", gsub("^.{7}", "", var_chr))
  if (any(vars == "veg_parm_ID")) {vars = vars[vars != "veg_parm_ID"]}
  
  value_match = regmatches(inlines, gregexpr("\\$1=.*,\\$2", inlines))
  value_chr = unlist(value_match[value_match != "character(0)"])
  values = gsub(".{3}$","", gsub("^.{3}", "", value_chr))
  
  return(data.frame(vars, values))
}
# ----------

shrub = get_vars("awk_spinup/init_shrubs.awk")
conifer = get_vars("awk_spinup/init_conifer.awk")
soils31 = get_vars("awk_spinup/init_soil31.awk")


# move files
file.copy("Preprocessing/Preprocess_Out/BC_conifer.world", to = "worldfiles/no_spin/BC_conifer.world", overwrite = TRUE)
file.copy("Preprocessing/Preprocess_Out/BC_shrub.world", to = "worldfiles/no_spin/BC_shrub.world", overwrite = TRUE)
file.copy("Preprocessing/Preprocess_Out/BC_shrub_conifer.world", to = "worldfiles/no_spin/BC_shrub_conifer.world", overwrite = TRUE)

file.copy("Preprocessing/Preprocess_Out/BC_conifer.flow", to = "flowtables/BC_conifer.flow", overwrite = TRUE)
file.copy("Preprocessing/Preprocess_Out/BC_shrub.flow", to = "flowtables/BC_shrub.flow", overwrite = TRUE)
file.copy("Preprocessing/Preprocess_Out/BC_shrub_conifer.flow", to = "flowtables/BC_shrub_conifer.flow", overwrite = TRUE)


worldfiles_in = c("worldfiles/no_spin/BC_shrub.world", "worldfiles/no_spin/BC_shrub_conifer.world", "worldfiles/no_spin/BC_conifer.world")

#shrub
update_world("worldfiles/no_spin/BC_shrub.world", "worldfiles/spin_soil_CN/BC_shrub.world", vars = shrub[,1], values = shrub[,2], overwrite = TRUE)
update_world("worldfiles/spin_soil_CN/BC_shrub.world","worldfiles/spin_soil_CN/BC_shrub.world", vars = soils31[,1], values = soils31[,2], overwrite = TRUE)

#conifer
update_world("worldfiles/no_spin/BC_conifer.world", "worldfiles/BC_conifer.world", vars = shrub[,1], values = shrub[,2], veg_parm_ID = 50, overwrite = TRUE)
update_world("worldfiles/spin_soil_CN/BC_conifer.world", "worldfiles/spin_soil_CN/BC_conifer.world", vars = conifer[,1], values = conifer[,2], veg_parm_ID = 7, overwrite = TRUE)
update_world("worldfiles/spin_soil_CN/BC_conifer.world","worldfiles/spin_soil_CN/BC_conifer.world", vars = soils31[,1], values = soils31[,2], overwrite = TRUE)

# conifer and shrub
update_world("worldfiles/no_spin/BC_shrub_conifer.world", "worldfiles/spin_soil_CN/BC_shrub_conifer.world", vars = shrub[,1], values = shrub[,2], veg_parm_ID = 50, overwrite = TRUE)
update_world("worldfiles/spin_soil_CN/BC_shrub_conifer.world", "worldfiles/spin_soil_CN/BC_shrub_conifer.world", vars = conifer[,1], values = conifer[,2], veg_parm_ID = 7, overwrite = TRUE)
update_world("worldfiles/spin_soil_CN/BC_shrub_conifer.world","worldfiles/spin_soil_CN/BC_shrub_conifer.world", vars = soils31[,1], values = soils31[,2], overwrite = TRUE)

worldfiles_in = file.path("worldfiles/spin_soil_CN/", list.files("worldfiles/spin_soil_CN", pattern = "BC_"))

# veg parm 50 = shrub cover fraction set to 0.6, gap fraction set to 0.0
# veg parm 7 = conifer cover fraction set to 0.6. gap fraction set to 0.1
for (f in worldfiles_in) {
  update_world(f, f, vars = c("cover_fraction", "gap_fraction"), values = c(0.6, 0), veg_parm_ID = 50, overwrite = TRUE)
  if (grepl("conifer", f) ) {
    update_world(f, f, vars = c("cover_fraction", "gap_fraction"), values = c(0.6, 0.1), veg_parm_ID = 7, overwrite = TRUE)
  }
}


# ----- Worldfile Changes -----
makeworlds = FALSE

if (makeworlds) {
  
  worldfiles = NULL
  # world files initialized with janets spun up values
  worldfiles_in = c("worldfiles/spin_soil_CN/BC_shrub.world", "worldfiles/spin_soil_CN/BC_shrub_conifer.world", "worldfiles/spin_soil_CN/BC_conifer.world")
  elevations = c(1000, 1500, 3000)
  aspect = c(90, 270)
  
  for (j in 1:length(worldfiles_in)) {
    out_name = sub(".world", "_N.world",worldfiles_in[j])
    worldfiles = c(worldfiles, out_name)
    update_world(worldfile = worldfiles_in[j], out_file = out_name, vars = c("z", "aspect"), values = c(elevations[j], 90), overwrite = TRUE)
    
    out_name = sub(".world", "_S.world",worldfiles_in[j])
    worldfiles = c(worldfiles, out_name)
    update_world(worldfiles_in[j], out_file = out_name, vars = c("z", "aspect"), values = c(elevations[j], 270), overwrite = TRUE)
  }
  
  for (i in 1:length(worldfiles)) { # to fix the wrong horizons
    update_world(worldfiles[i], out_file = worldfiles[i], vars = c("e_horizon", "w_horizon"), values = c(0.08715574016, 0.2079116808), overwrite = TRUE)
  }
  
  #change cover fraction for all vegs
  for (i in 1:length(worldfiles)) { # to fix the wrong horizons
    update_world(worldfiles[i], out_file = worldfiles[i], vars = c("cover_fraction"), values = c(0.6), overwrite = TRUE)
  }
  
  
  #change gap fraction for conifer
  change_gap = c(
    "worldfiles/BC_shrub_conifer_N.world",
    "worldfiles/BC_shrub_conifer_S.world",
    "worldfiles/BC_conifer_N.world",
    "worldfiles/BC_conifer_S.world"
  )
  for (i in seq_along(change_gap)) {
    update_world(worldfile = change_gap[i], out_file = change_gap[i], vars = c("gap_fraction"), values = c(0.1), veg_parm_ID = "7", overwrite = TRUE)
  }
  
}




