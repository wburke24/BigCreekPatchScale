setwdhere()

inst = FALSE
if (inst) {
  devtools::install("../../../rhessys/RHESSysIOinR")
  devtools::install_github("RHESSys/RHESSysIOinR")
  devtools::install_github("RHESSys/RHESSysPreprocessing")
}

library(RHESSysIOinR)
library(RHESSysPreprocessing)

# any will do 
worldfiles_in = c("worldfiles/spin_soil_CN/BC_conifer_N.world", "worldfiles/spin_soil_CN/BC_shrub_conifer_N.world", "worldfiles/spin_soil_CN/BC_shrub_N.world")

scenarios = NULL
load("scenarios.rdata")

# ----- generate redefines for treatments/thinning -----
# basename for output redefs
redef_worldfiles = file.path("worldfiles/redefine/", gsub("_N", "", basename(worldfiles_in)))

for (w in seq_along(worldfiles_in)) {
  
  # if for 8 asppatches
  if ( grepl("shrub_conifer", worldfiles_in[w]) ) {
    
    # 101 102 103 104 105 106 107 108 109
    # pct_family_area	value 0.05 | 0.05 | 0.125 | 0.125 | .2 | .2 | 0.05 | .075 | 0.125
    # _canopy_strata 2 | 1 | 2 | 1 | 2 | 1 | 2 | 2 | 1 
    # veg_parm_ID		dvalue 7  50 | 50 | 7  50 | 50 | 7  50 | 50 | 7 50 | 7 50 | 50 
    
    # remove pct understory + lagged remove all litter - for redefine_world_thin_X remain
    build_redefine(worldfiles_in[w], sub(".world", "_under_0.1.world", redef_worldfiles[w]), std_thin = "1", veg_parm_ID = "50", patchID = c("101", "102"))
    
    
    
    build_redefine(worldfiles_in[w], sub(".world", "_under_0.25.world" ,redef_worldfiles[w]),  std_thin = "1", veg_parm_ID = "50", patchID = c("103", "104"))
    build_redefine(worldfiles_in[w], sub(".world", "_under_0.4.world" ,redef_worldfiles[w]),  std_thin = "1", veg_parm_ID = "50", patchID = c("105", "106"))
    
    # litter to 0 <patch level> and cwd to 0 <strata> - for redefine_world - combine w above for first 3 thin options and w remove understory for prescribed fire
    litter_vars = c("litter_cs.litr1c", "litter_ns.litr1n",	"litter_cs.litr2c", "litter_cs.litr3c", "litter_cs.litr4c", "cs.cwdc", "ns.cwdn")
    build_redefine(worldfiles_in[w], sub(".world", "_litter_0.world", redef_worldfiles[w]), vars = litter_vars, values = 0)
    
    # selection thinning (remove pct overstory) + carbon to litter/carbon remains - redefine_world_thin_X (remain/harvest)
    build_redefine(worldfiles_in[w], sub(".world", "_over_0.1.world" ,redef_worldfiles[w]),  std_thin = "1", veg_parm_ID = "7", patchID = c("101", "107"))
    build_redefine(worldfiles_in[w], sub(".world", "_over_0.25.world" ,redef_worldfiles[w]),  std_thin = "1", veg_parm_ID = "7", patchID = c("103", "107", "108"))
    build_redefine(worldfiles_in[w], sub(".world", "_over_0.4.world" ,redef_worldfiles[w]),  std_thin = "1", veg_parm_ID = "7", patchID = c("105", "103", "108"))
    
    # thin understory - all patches, by half, combine w litter to 0, for prescribed fire
    #build_redefine(worldfilesB[w], sub(".world", "_under_0.5_allpatch.world" ,redef_worldfilesB[w]),  std_thin = ".5", veg_parm_ID = "50")
    
  } else {
    # pct_family_area value .10 | .25 | .40 | .25
    # _canopy_strata 1 | 1 | 1 | 1
    # veg_parm_ID		dvalue 50 | 50 | 50 | 50 
    
    # remove pct understory + lagged remove all litter - for redefine_world_thin_X remain
    # done via removal of shrubs from associated coverage area (aspatial patch), eg patch w 10% coverage, etc
    # std_thin/value = fraction to be removed
    build_redefine(worldfiles_in[w], sub(".world", "_under_0.1.world" ,redef_worldfiles[w]),  std_thin = "1", veg_parm_ID = "50", patchID = "101")
    build_redefine(worldfiles_in[w], sub(".world", "_under_0.25.world" ,redef_worldfiles[w]),  std_thin = "1", veg_parm_ID = "50", patchID = "102")
    build_redefine(worldfiles_in[w], sub(".world", "_under_0.4.world" ,redef_worldfiles[w]),  std_thin = "1", veg_parm_ID = "50", patchID = "103")
    
    # litter to 0 <patch level> and cwd to 0 <strata> - for redefine_world - combine w above for first 3 thin options and w remove understory for prescribed fire
    litter_vars = c("litter_cs.litr1c", "litter_ns.litr1n",	"litter_cs.litr2c", "litter_cs.litr3c", "litter_cs.litr4c", "cs.cwdc", "ns.cwdn")
    build_redefine(worldfiles_in[w], sub(".world", "_litter_0.world", redef_worldfiles[w]), vars = litter_vars, values = 0)
    
    # selection thinning (remove pct overstory) + carbon to litter/carbon remains - redefine_world_thin_X (remain/harvest)
    # this doesn't make sense for shrub only options and doesn't do anything rn - uneeded runs are removed in the scenarios table
    build_redefine(worldfiles_in[w], sub(".world", "_over_0.1.world" ,redef_worldfiles[w]),  std_thin = "1", veg_parm_ID = "7", patchID = "101")
    build_redefine(worldfiles_in[w], sub(".world", "_over_0.25.world" ,redef_worldfiles[w]),  std_thin = "1", veg_parm_ID = "7", patchID = "102")
    build_redefine(worldfiles_in[w], sub(".world", "_over_0.4.world" ,redef_worldfiles[w]),  std_thin = "1", veg_parm_ID = "7", patchID = "103")
    
    # thin understory - all patches, by half, combine w litter to 0, for prescribed fire
    #build_redefine(worldfilesA[w], sub(".world", "_under_0.5_allpatch.world" ,redef_worldfilesA[w]),  std_thin = ".5", veg_parm_ID = "50")
  }
  
  
  
} # end loop through worldfiles


# remove shrub + overstory thinning to be sure they won't be used in error
file.remove(c("worldfiles/redefine/BC_shrub_over_0.1.world", "worldfiles/redefine//BC_shrub_over_0.25.world", "worldfiles/redefine//BC_shrub_over_0.4.world"))




