# Big Creek workflow

setwdhere()

# Make the scenarios
source("BCp_1_make_scenarios.R")

# make/get the correct worldfiles - modify as needed
source("BCp_2_make_worldfiles.R")

# spin up
source("BCp_3_spinup.R")

# make redefines - only needed to re run if name changes have happened
source("BCp_4_make_redefines.R")

# run scenarios
source("BCp_5_run.R")

# collect output - collect from patches to patch families
source("BCp_6_collect_output.R")

# processes output / create plots
source("BCp_7_output.R")