# check soil water storage

# TEST/CHECK STORAGE - see img in gdrive
fc * rootdepth/ * (1 - sat-to-gw (GW1 but in def file version))
fc = 0.5 * (pa(psi air entry) /3.4) ^ po(pore size)
fc = 0.5 *  (input_def_list[[28]][[3]] / 3.4) ^ input_def_list[[26]][[3]]
storage = fc * input_def_list[[25]][[3]] * (1 - input_def_list[[27]][[3]])