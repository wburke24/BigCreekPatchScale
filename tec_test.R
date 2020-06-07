if (scenarios$treatments[i] != "no treatment") {
  treat = unlist(strsplit(scenarios$treatments[i], "\\s"))
  worldfiles_base = substr(scenarios$worldfile[i], 12, nchar(scenarios$worldfile[i]) - 6)
  redef = file.path("worldfiles/redefine", list.files("worldfiles/redefine/", paste0(worldfiles_base, "_", treat[1])))
  
  if (startsWith(scenarios$treatments[i], "under")) {
    # this one has 2 treatments - thin from below and prescribed fire
    tec_type = "redefine_world_thin_remain"
    tec_treat = tec_repeat(start, end, repeat_interval[scenarios$interval[i]], "year", tec_type, scenarios$worldfile[i], redef, overwrite = TRUE)
    redef = file.path("worldfiles/redefine", list.files("worldfiles/redefine/", paste0(worldfiles_base, "_litter_0")))
    inc1mo = seq(as.Date(start, "%Y %m %d %H"), by = "month", length = 2)[2]
    start1mo = paste0(gsub("-"," ", inc1mo), " 10")
    tec_type = "redefine_world"
    tec_treat2 = tec_repeat(start1mo, end, repeat_interval[scenarios$interval[i]], "year", tec_type, scenarios$worldfile[i], redef, overwrite = TRUE)
    tec_treat = rbind(tec_treat, tec_treat2)
    #tec_type = "redefine_world"
    #tec_treat = tec_repeat(start, end, repeat_interval[scenarios$interval[i]], "year", tec_type, scenarios$worldfile[i], redef)
  } else if (startsWith(scenarios$treatments[i], "over") & endsWith(scenarios$treatments[i], "remains")) {
    tec_type = "redefine_world_thin_remain"
    tec_treat = tec_repeat(start, end, repeat_interval[scenarios$interval[i]], "year", tec_type, scenarios$worldfile[i], redef, overwrite = TRUE)
  } else if (startsWith(scenarios$treatments[i], "over") & endsWith(scenarios$treatments[i], "removed")) {
    tec_type = "redefine_world_thin_harvest"
    tec_treat = tec_repeat(start, end, repeat_interval[scenarios$interval[i]], "year", tec_type, scenarios$worldfile[i], redef, overwrite = TRUE)
  } else if (startsWith(scenarios$treatments[i], "litter")) {
    tec_type = "redefine_world"
    tec_treat = tec_repeat(start, end, repeat_interval[scenarios$interval[i]], "year", tec_type, scenarios$worldfile[i], redef, overwrite = TRUE)
  }
  input_tec_data = rbind(input_tec_data, tec_treat)
  input_tec_data = input_tec_data[with(input_tec_data, order(year, month, day, hour)), ]
  # maybe test if there's an overlap
}




test = match(unique(scenarios$treatments),scenarios$treatments)
scenarios[test,]

runct = 0

for (i in test) {
  
  runct = runct + 1
  cat("\n ========== Scenario ",runct," ========== \n\n")
  
  # ----- input_rhessys vars -----
  input_rhessys$start_date <- scenarios$start_date[i]
  input_rhessys$end_date <- scenarios$end_date[i]
  input_rhessys$world_file = scenarios$worldfile[i]
  input_rhessys$flow_file = sub("worldfiles", "flowtables", sub("(_N|_S)\\.world",".flow", scenarios$worldfile[i]))
  
  # ----- veg, soil and sharing parameters -----
  input_def_list = c(list(list(input_hdr_list$landuse_def[1], "sh_l", scenarios$sharing[i])),soilslist[[scenarios$soils[i]]])
  #input_def_list = c(list(list(input_hdr_list$landuse_def[1], "sh_l", 0.25)),soilslist[[scenarios$soils[i]]])
  #list(input_hdr_list$stratum_def[1], "epc.height_to_stem_coef", 9.4),
  #list(input_hdr_list$stratum_def[1], "epc.height_to_stem_coef", 18.39)), 
  
  # ----- climate/CO2 -----
  input_clim_base_list[[1]]$daily$c1[1] = scenarios$climate[i]
  if (scenarios$climate[i] == "clim/upGGmod") {
    input_def_list = c(input_def_list, list(list(input_hdr_list$zone_def[1], "atm_CO2", 450)))
  }
  
  # ----- treatments - tec events -----
  input_tec_data = input_tec(start = scenarios$start_date[i], end = scenarios$end_date[i], output_state = FALSE)
  
  # thinning type - NA is no thinning (could still have prescribed fire)
  if (!is.na(scenarios$treat_type[i])) {
    redef_world = paste0(substr(scenarios$worldfile[i], 12, nchar(scenarios$worldfile[i]) - 6), "_", 
                         scenarios$treat_type[i], "_", scenarios$treat_intensity[i], ".world")
    redef_path = file.path("worldfiles/redefine", redef_world)
    
    tec_treat = tec_repeat(scenarios$start_date[i], scenarios$end_date[i], scenarios$interval[i], "year", 
                           paste0("redefine_world_", scenarios$tec_type[i]), scenarios$worldfile[i], redef_path, overwrite = TRUE)
    
    input_tec_data = rbind(input_tec_data, tec_treat)
  }
  if (scenarios$presc_fire[i]) {
    redef_world = paste0(substr(scenarios$worldfile[i], 12, nchar(scenarios$worldfile[i]) - 6), "_litter_0.world")
    redef_path = file.path("worldfiles/redefine", redef_world)
    start_1mo = paste0(gsub("-"," ", seq(as.Date(scenarios$start_date[i], "%Y %m %d %H"), by = "month", length = 2)[2]), " 10")
    
    tec_presc = tec_repeat(start_1mo, scenarios$end_date[i], scenarios$interval[i], "year", "redefine_world", 
                           scenarios$worldfile[i], redef_path, overwrite = TRUE)
    
    input_tec_data = rbind(input_tec_data, tec_presc)
  }
  input_tec_data = input_tec_data[with(input_tec_data, order(year, month, day, hour)), ]

  print(input_tec_data)
    
}




